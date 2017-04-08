import numpy as np
from astropy import units as u
from astropy import constants as c
import pandas as pd
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import  RectBivariateSpline
from scipy.stats import binom
from scipy.interpolate import interp1d
from matplotlib.pylab import * 
from scipy.optimize import newton

TDUR_EARTH_SUN_HRS = (
    ((4 * c.R_sun**3 * 1.0*u.yr / np.pi / c.G / (1.0*c.M_sun))**(1.0/3.0)).to(u.hr)).value
DEPTH_EARTH_SUN = ((c.R_earth / c.R_sun)**2).cgs.value


def _sbins_to_sbins_id(sbins):
    """Convenience function for converting bin names to bins_id names"""
    return ['{}_bin_id'.format(sbin) for sbin in sbins]

class Grid(object):
    """
    Grid

    This object creates grids used to compute occurrence and completeness
    """

    def __init__(self, bins_dict, spacing_dict):
        """Lay down grid for computing occurrence
        
        Args:
            bins_dict (dict): Dictionary of lists defining bins. e.g.:
                {'per': [5,10,20], 'prad': [1,2,4], 'smet': [-0.4,0.0,0.4]}
            spacing_dict (dict): Specify linear or log spacing e.g.:
                {'per': 'log', 'prad': 'log', 'smet': 'linear'}
        
        """

        assert bins_dict.keys()==spacing_dict.keys()
        self.bins = bins_dict
        self.spacing = spacing_dict

        # For every bin, record the lower, upper, and central value
        self.bins1 = {}
        self.bins2 = {}
        self.binsc = {}
        self.nbins = {}

        for key in self.bins.keys():
            self._add_bin_edges(key)
            self.nbins[key] = len(self.bins1[key])


        ds = xr.Dataset(coords=self.binsc)
        for sbin in self.bins.keys():
            ds.coords[sbin] = self.binsc[sbin]

        for key in self.bins.keys():
            coord = {key: self.binsc[key]}
            ds[key+'c'] = xr.DataArray(self.binsc[key],coord)
            ds[key+'1'] = xr.DataArray(self.bins1[key],coord)
            ds[key+'2'] = xr.DataArray(self.bins2[key],coord)

        # This line makes all data variables 2D
        ds = ds.to_dataframe().to_xarray()
        self.ds = ds

    def _add_bin_edges(self, key):
        """Computes bin edges and centers for the specified bins"""
        bins = self.bins[key]
        spacing = self.spacing[key]

        bins1 = []
        bins2 = []
        binsc = []
        for i in range(len(bins)-1):
            bin1 = bins[i]
            bin2 = bins[i+1]
            if spacing=='linear':
                binc = 0.5 * (bin1 + bin2)
            elif spacing=='log':
                binc = np.sqrt(bin1 * bin2)
            else:
                False, "Invalid spacing"
            bins1+=[bin1]
            bins2+=[bin2]
            binsc+=[binc]
        
        self.bins1[key] = bins1
        self.bins2[key] = bins2
        self.binsc[key] = binsc

class Completeness(object):
    """
    Class for that can compute completeness using the noise properties
    of a representive ensemble of stars.  
    """

    def __init__(self, stars, grid, prob_det_mes_name, impact_transit):
        """
        Args:
            stars (pandas.DataFrame): Sample of stars from which planets 
                are detected. Must be as close as possible to be sample of 
                of stars used in the planet search. Must contain the following
                keys.
                - logcdpp3: three hour CDPP
                - logcdpp6: six
                - logcdpp12: twelve
                - tobserved: days in which target was observed
                - smass: Stellar mass (solar masses) 
                - srad: Stellar radius (solar-radii) 
            grid (Grid): Grid object that contains boundaries of bins.

        """
        self.stars = stars

        # Define the regular grid for interpolation
        x0 = np.array(stars.index) # (nstar)
        x1 = np.log10([3,6,12]) # (3)
        points  = (x0, x1) 
        cols = 'logcdpp3 logcdpp6 logcdpp12'.split()
        values = stars[cols] # (nstar, 3)
        values = np.array(values)
        logcdpp_interp = RegularGridInterpolator(
            points, values, bounds_error=False, fill_value=None
        )        

        self.logcdpp_interp = logcdpp_interp
        self.x0 = x0 
        self.x1 = x1
        self.values = values
        self.nstars = len(stars)
        self.grid = grid
        self.prob_det_mes_name = prob_det_mes_name
        self.impact_transit = impact_transit

    def mes(self, per, prad):
        """
        Calculate Multiple Event Statistic

        For a planet of a given period and size, calculate the
        expected MES for every star in the sample.

        Args:
            per (float) : orbital period of planet days
            prad (float): size of planet in Earth-radii 

        Returns:
            array: MES

        """
        depth = self._depth(prad)
        tdur = self._tdur(per)
        cdpp = self._cdpp(tdur) * 1e-6
        num_transits = self._num_transits(per)
        _mes = depth / cdpp * np.sqrt(num_transits)
        return _mes

    def mes_scaled(self, per, prad):
        """
        Calculate Multiple Event Statistic and apply scaling
        """
        return self.mesfac * self.mes(per, prad)

    def _tdur(self, per):
        """
        Compute duration for a putative transit of a given orbital period
        
        Args:
            per (float): orbital period

        Returns:
            pandas.Series: transit duration for each star in the sample.
        """
        per_yrs = per / 365.25
        tdur = ( 
            TDUR_EARTH_SUN_HRS * self.stars['srad'] * 
            (per_yrs / self.stars['smass'] )**(1.0/3.0)
        )
        return tdur

    def _smax(self, per):
        """
        Compute semi-major axis putative transit of a given orbital period
        
        Args:
            per (float): orbital period

        Returns:
            pandas.Series: transit duration for each star in the sample.
        """

        
        per_yrs = per / 365.25
        smax = self.stars['smass']**(1.0/3.0) * per_yrs**(2.0/3.0)
        return smax

    def _depth(self, prad):
        """
        Compute transit depth for a putative planet of a given size.
        
        Args: 
            prad (float): planet radius (Earth-radii) 

        Returns:
            pd.Series: the depth of a `prad` planet around each star

        """
        _depth = (prad / self.stars.srad)**2 * DEPTH_EARTH_SUN
        return _depth

    def _num_transits(self, per):
        """
        Compute number of transits for a planet having a given orbital
        period, factoring in duty cycle.

        Args:
            per (float): orbital period (days)
          
        Returns:
            pd.Series: the number of transits for each star in the sample
        """

        return self.stars.tobs / per

    def _cdpp(self, tdur):
        """
        Calculate noise (CDPP) over a specified duration
        
        Args:
            tdur (pandas.Series or array): the transit duration for each star

        Retruns:
            pd.Series: The CDPP for each star

        """
        tdur = np.array(tdur)
        logtdur = np.log10(tdur)

        x0i = self.x0
        x1i = logtdur
        pointsi = np.vstack([x0i,x1i]).T
        logcdpp = self.logcdpp_interp(pointsi)
        cdpp = pd.Series(index=self.stars.index,data=10**logcdpp)
        return cdpp

    def prob_det_mes(self, mes):
        """Recovery rate vs. multiple event statistic

        Args:
            mes (array): Multiple event statistic
            
        """
        out = np.zeros_like(mes) 
        mes = np.array(mes) 
        if self.prob_det_mes_name.count('step')==1:
            min_mes = self.prob_det_mes_name.split('-')[1]
            min_mes = float(min_mes)
            out[mes > min_mes] = 1.0

        return out

    def prob_det(self, per, prad, method='direct'):
        """Probability that a planet would be detectable

        Probability that transiting planet with orbital period `per`
        and planet radius `prad` would be detectable around a
        randomly-drawn star from the stellar sample.

        Args:
            per (float): orbital period
            prad (float): planet size
            method (str): One of the following:
                - direct: compute prob of detectability by counting stars 
                  where the MES is sufficient to detect.
                - interp: use the interpolator

        Returns:
            float: fraction of stars in sample where we could have detected 
                planet

        """        
        
        if method=='direct':
            mes = self.mes_scaled(per, prad)
            _prob_det = 1.0*self.prob_det_mes(mes).sum() / self.nstars 
        elif method=='interp': 
            #_prob_det = self.prob_det_interp((per, prad))
            _prob_det = self.prob_det_interp(per, prad, grid=False)
        else:
            assert False,'Invalid method'

        return _prob_det

    def prob_tr(self, per):
        """
        Probability that a given planet would transit
        """
        a = self._smax(per)
        srad = (np.array(self.stars['srad']) * u.R_sun).to(u.AU).value
        _prob_tr = srad * self.impact_transit / a
        return _prob_tr.mean()

    def compute_mes_factor(self, plnt_mes):
        """
        Using a simple SNR calculation, we tend to over estimate the MES
        that the pipeline would have found.

        Args:
            plnt_mes (pd.DataFrame): Planets used to calibrate MES. Must have 
                the fllowing keys:
                - per
                - prad
                - id_kic
                - id_koid
             
        """
        plnt_mes = plnt_mes.copy()
        self.plnt_mes = plnt_mes
        self.plnt_mes['mes_pipeline'] = self.plnt_mes.mes
        self._compute_mes_factor_loop(False)

        logmes_pipeline = np.log10(self.plnt_mes.mes_pipeline)
        logmes_formula = np.log10(self.plnt_mes.mes_formula)
        logmes_diff_med = np.nanmedian(logmes_formula - logmes_pipeline)
        logmes_diff_std = np.nanstd(logmes_formula - logmes_pipeline)
        print "med(log(mes_formula/mes_pipeline)) {:.2f} (dex)".format(
            logmes_diff_med
        )
        self.logmes_diff_std = logmes_diff_med
        self.logmes_diff_med = logmes_diff_std 
        self.mesfac = 10**(-1.0 * logmes_diff_med)
        self._compute_mes_factor_loop(True)
 
    def _compute_mes_factor_loop(self, scaled):
        i = 0
        for id_koicand, row in self.plnt_mes.iterrows():
            try:
                if scaled:
                    fmes = self.mes_scaled
                    key = 'mes_formula_scaled' 
                else:
                    fmes = self.mes
                    key = 'mes_formula' 

                mes_formula = fmes(row.per,row.prad).ix[row.id_kic]
                mes_pipeline = row.mes
                self.plnt_mes.ix[id_koicand,key] = mes_formula

                if i< 20:
                    s =  "{} {:.2f} {:.2f}".format(
                        id_koicand, mes_pipeline, mes_formula
                    )
                    print s
            except KeyError:
                pass
        
            i+=1

    def compute_grid_prob_det(self,verbose=0):
        """Compute a grid of detection probabilities"""
        def rowfunc(row):
            return self.prob_det(row.perc, row.pradc)

        def callback(row):
            s ="{perc:.3f} {pradc:.3f} {temp:.3f}".format(**row)
            if verbose>0:
                print s

        self.grid.ds['prob_det'] = self._grid_loop(rowfunc, callback)

    def compute_grid_prob_tr(self,verbose=0):
        """Compute a grid of transit probabilities"""
        def rowfunc(row):
            return self.prob_tr(row.perc)

        def callback(row):
            s ="{perc:.3f} {pradc:.3f} {temp:.3f}".format(**row)
            if verbose>0:
                print s

        self.grid.ds['prob_tr'] = self._grid_loop(rowfunc, callback)

    def _grid_loop(self, rowfunc, callback):
        df = self.grid.ds.to_dataframe()
        i = 0
        for idx, row in df.iterrows():
            df.ix[idx,'temp'] = rowfunc(row)
            if i%100==0:
                callback(df.ix[idx])
            i+=1
        return df['temp'].to_xarray()

    def init_prob_det_interp(self):
        """
        
        """
        # Define the regular grid for interpolation
        x0 = np.array(self.grid.ds.per)
        x1 = np.array(self.grid.ds.prad)
        points  = (x0, x1) 
        values = self.grid.ds['prob_det'].transpose('per','prad')
        #self.prob_det_interp = RegularGridInterpolator(
        #    points, values, bounds_error=False
        #)        
        self.prob_det_interp = RectBivariateSpline(
            x0, x1, values
        )        

class Occurrence(object):
    """
    Occurrence Object

    This object allows one to compute planet occurrence in bins of
    orbital period, planet radius, and stellar metallicity

    Args:
        plnt (pandas.DataFrame): One of each row per planet
            - srad: stellar radius (Solar-radii)
            - per: planet orbital period (days)
            - prad: planet radius (Earth-radii) 
            - smax: planet semi-major axis (AU)
            - mes: planet multiple event statistic
            - impact: planet impact parameter
        nstars (float): number of stars from which our sample was drawn
        comp (Completeness): completeness object
    
    """
    plnt_columns = ['srad','per','prad','smax','impact','mes']

    def __init__(self, plnt, nstars, comp, grid):
        for col in self.plnt_columns:
            assert list(plnt).count(col)==1, \
                "DataFrame must contain {}".format(col)

        idxdrop = plnt[plnt.impact > comp.impact_transit].index
        ndrop = len(idxdrop)
        print "dropping {} planets, impact > {}".format(
            ndrop, comp.impact_transit
        )

        plnt = plnt.drop(idxdrop)
        if comp.prob_det_mes_name.count('step')==1:
            min_mes = comp.prob_det_mes_name.split('-')[1]
            min_mes = float(min_mes)
            idxdrop = plnt[plnt.mes < min_mes].index
            ndrop = len(idxdrop)
            print "dropping {} planets, mes < {}".format(ndrop,min_mes)

        plnt = plnt.drop(idxdrop)

        self.plnt = plnt
        self.nstars = nstars
        self.comp = comp

    def set_grid_ntrials(self, grid):
        """
        Compute number of effective trials in little boxes in
        period-radius space

        Args:
            grid (Grid object): grid small compared to variability in 
               completeness
        """

        # Convert to dataframe so we can loop over
        df = grid.ds.to_dataframe()
        for i, row in df.iterrows():
            prob_det = self.comp.prob_det(row.perc,row.pradc,method='direct')
            prob_tr = self.comp.prob_tr(row.perc)
            df.ix[i,'prob_det'] = prob_det
            df.ix[i,'prob_tr'] = prob_tr
            
        df['ntrials'] = self.nstars * df['prob_det'] * df['prob_tr']
        grid.ds = df.to_xarray()
        self.grid_ntrials = grid




    def compute_occurrence(self):
        """Compute occurrence over grid

        For a planet at the center of the bin what fraction of the
        stars in the sample, could we detect?

        Args:
            sbins (list): list of strings specifying the bins to compute 
                occurrence over

        Returns:
            xarray: data with the following quantities computed in each bin
                - plnt_transit_sum: number of transiting planets
                - plnt_total_sum: number planets (including non-transiting)
        """
        plnt = self.plnt
        a = np.array(plnt['smax']) * u.AU
        srad = np.array(plnt['srad']) * c.R_sun

        # Compute detectability for each transit
        plnt['prob_tr'] = (srad * self.comp.impact_transit / a).cgs.value
        plnt['prob_det'] = self.comp.prob_det(
            plnt['per'], plnt['prad'], method='interp'
        )
        plnt['prob_trdet'] = plnt['prob_tr'] * plnt['prob_det']
        for suffix in 'prob_tr prob_det prob_trdet'.split():
            plnt['inv_'+suffix] = plnt[suffix]**-1.0

        # Compute agg quantities in bins
        labels = []
        for sbin in self.grid.bins.keys():
            labels+=[pd.cut(self.plnt[sbin],self.grid.bins[sbin])]

        g = plnt.groupby(labels,as_index=True)
        df = g[['id_starname']].count()
        
        # Transiting and detected planets per bin
        df['plnt_trdet_sum'] = g[['inv_prob_tr']].count()

        # Transiting planets per bin (accounting for missed transiting planets)
        df['plnt_tr_sum'] = g[['inv_prob_det']].sum()

        # Planets per bin (accounting for missed and non-transiting planets)
        df['plnt_sum'] = g[['inv_prob_trdet']].sum()
        df['plnt_occur'] = df['plnt_sum'] / self.nstars
        ds = df.to_xarray()
        for sbin in self.grid.bins.keys():
            ds.coords[sbin] = self.grid.binsc[sbin]

        ds = xr.merge([self.grid.ds,ds])
        return ds 

    def compute_occurrence_errors(self):
        """
        Compute occurrence with errors
        
        Look for planets in each bin. If they exist, compute 1000 samples of the pdf
        """
        # check that a/Rstar does not change by more than 10% over a bin
        # check that prob_detectr does not change by more than 10% over a bin

        # Compute agg quantities in bins
        prob_tr_const_thresh = 0.2 # 
        prob_det_const_thresh = 0.2 # 

        labels = []
        for sbin in self.grid.bins.keys():
            labels+=[pd.cut(self.plnt[sbin],self.grid.bins[sbin])]

        g = self.plnt.groupby(labels,as_index=True)
        df = g[['id_starname']].count()
        ds = df.to_xarray()
        for sbin in self.grid.bins.keys():
            ds.coords[sbin] = self.grid.binsc[sbin]

        ds = xr.merge([self.grid.ds,ds])
        df = ds.to_dataframe()
        df['is_prob_tr_const'] = 0
        df['is_prob_det_const'] = 0
        df['prob_tr_lo'] = 0
        df['prob_tr_hi'] = 0

        for i, row in df.iterrows():
            # Difference in transit prob
            prob_tr_hi = self.comp.prob_tr(row.per1)
            prob_tr_lo = self.comp.prob_tr(row.per2)
            df.ix[i,'prob_tr_lo'] = prob_tr_lo
            df.ix[i,'prob_tr_hi'] = prob_tr_hi
            
            # Difference in completeness
            prob_det_hi = self.comp.prob_det(row.per1,row.prad2,method='direct')
            prob_det_lo = self.comp.prob_det(row.per2,row.prad1,method='direct')
            df.ix[i,'prob_det_lo'] = prob_det_lo
            df.ix[i,'prob_det_hi'] = prob_det_hi

            ratio = prob_tr_hi/prob_tr_lo
            if np.abs(ratio-1) < prob_tr_const_thresh:
                df.ix[i,'is_prob_tr_const'] = 1

            ratio = prob_det_hi/prob_det_lo
            # print row.per1, row.prad2, row.per2, row.prad1, prob_det_hi, prob_det_lo
            if np.abs(ratio-1) < prob_det_const_thresh:
                df.ix[i,'is_prob_det_const'] = 1

        pdfL = []
        rateL = []
        df['rate_ul'] = -1
        df['prob_det'] = -1
        df['rate_prob_tr'] = -1
        for i, row in df.iterrows():
            cut = self.plnt[
                self.plnt.prad.between(row.prad1,row.prad2) & 
                self.plnt.per.between(row.per1,row.per2)
            ]
            nplanets = len(cut)
            nstars = round(self.nstars)

            prob_tr = self.comp.prob_tr(row.perc)
            prob_det = self.comp.prob_det(row.perc,row.pradc,method='interp')
            ntrials = np.round(self.nstars * prob_det)
            df.ix[i,'ntrials'] = ntrials
            df.ix[i,'nstars'] = nstars
            df.ix[i,'nplanets'] = nplanets
            df.ix[i,'prob_det'] = prob_det
            df.ix[i,'prob_tr'] = prob_tr

            if nplanets!=0:
                # Detection, calculate pdf of the rate of planets in box
                rate_tr_best = nplanets / ntrials
                rate_tr = np.linspace(0,10 * rate_tr_best,1000)
                pdf = [binom.pmf(nplanets, ntrials, r) for r in rate_tr]
                rate = rate_tr / prob_tr
                pdf /= np.sum(pdf)
            else:
                # Non-detection, calculate 95% upper limit on rate.
                obj = lambda rate : binom.pmf(0, ntrials, rate) - 0.05
                rate_tr = newton(obj,0)
                rate = rate_tr / prob_tr
                df.ix[i,'rate_ul'] = rate
                pdf = []
                rate = []

            pdfL.append(pdf)
            rateL.append(rate)

        df['pdf'] = pdfL
        df['rate'] = rateL
        return df

def combine_cells(df,plot_diag=False):
    """
    Combine cells

    Each cell contains either an upper limit, or a PDF of the rate. We
    either sum up all the rate PDFs (using convolution), or we take
    the maximum upper limit as the total upper limit.

    """

    # Figure out how much to interpolate
    idx =  df[df.nplanets>0].index

    out = {'rate_ul':None,'rate':None,'rate_err1':None,'rate_err2':None}
    if len(idx)==0:
        out['rate_ul'] = df['rate_ul'].max()
        return out

    rate2d = np.vstack(df.ix[idx,'rate'].values)
    pdf2d = np.vstack(df.ix[idx,'pdf'].values)
    rateistep = np.min(rate2d[:,1:] - rate2d[:,:-1])
    rateimax = np.max(rate2d)
    ratei = np.arange(0,rateimax,rateistep)

    # resample the rate onto a constant scale
    pdfi2d = []
    for i in range(pdf2d.shape[0]):
        rate = rate2d[i]
        pdf = pdf2d[i]
        interp = interp1d(rate, pdf, fill_value=0, bounds_error=False)
        pdfi = interp(ratei)
        pdfi2d.append(pdfi)

    pdfconv = reduce(np.convolve,pdfi2d)
    pdfconv /= np.sum(pdfconv)
    rateiconv = arange(len(pdfconv))*rateistep

    cdfconv = np.cumsum(pdfconv)
    ppf = interp1d(cdfconv,rateiconv)

    ratelo,ratemi,ratehi = ppf([.16,.500,0.84])
    rate_err1 = ratehi-ratemi
    rate_err2 = ratelo-ratemi

    rates = np.array([ratemi,rate_err1,rate_err2]) * 100
    print "{0:.3f} +/- {1:.3f} {2:.3f} %".format(*rates)

    if plot_diag:
        fig,axL = subplots(nrows=2,sharex=True)
        sca(axL[0])
        plot(rate2d.T,pdf2d.T,'k')
        plot(ratei,np.vstack(pdfi2d).T,'b')
        sca(axL[1])
        plot(rateiconv,pdfconv,'r')
        x = [ratemi]
        y = [1.1*max(pdfconv)]
        xerr = [[-rate_err2],[rate_err1]]
        errorbar(x, y, xerr=xerr, fmt='o')

    out['rate'] = ratemi
    out['rate_err1'] = rate_err1
    out['rate_err2'] = rate_err2
    return out

def combine_cells2(df,plot_diag=False):
    """
    Combine cells

    Each cell contains either an upper limit, or a PDF of the rate. We
    either sum up all the rate PDFs (using convolution), or we take
    the maximum upper limit as the total upper limit.

    """

    # Figure out how much to interpolate
    idx =  df[df.nplanets>0].index

    out = {'rate_ul':None,'rate':None,'rate_err1':None,'rate_err2':None}
    if len(idx)==0:
        out['rate_ul'] = df['rate_ul'].max()
        return out

    rate2d = np.vstack(df.ix[idx,'rate'].values)
    pdf2d = np.vstack(df.ix[idx,'pdf'].values)
    rateistep = np.min(rate2d[:,1:] - rate2d[:,:-1])
    rateimax = np.max(rate2d)
    ratei = np.arange(0,rateimax,rateistep)

    # resample the rate onto a constant scale
    pdfi2d = []
    for i in range(pdf2d.shape[0]):
        rate = rate2d[i]
        pdf = pdf2d[i]
        interp = interp1d(rate, pdf, fill_value=0, bounds_error=False)
        pdfi = interp(ratei)
        pdfi2d.append(pdfi)

    pdfconv = reduce(np.convolve,pdfi2d)
    pdfconv /= np.sum(pdfconv)
    rateiconv = arange(len(pdfconv))*rateistep

    cdfconv = np.cumsum(pdfconv)
    ppf = interp1d(cdfconv,rateiconv)

    ratelo,ratemi,ratehi = ppf([.16,.500,0.84])
    rate_err1 = ratehi-ratemi
    rate_err2 = ratelo-ratemi

    rates = np.array([ratemi,rate_err1,rate_err2]) * 100
    print "{0:.3f} +/- {1:.3f} {2:.3f} %".format(*rates)

    if plot_diag:
        fig,axL = subplots(nrows=2,sharex=True)
        sca(axL[0])
        plot(rate2d.T,pdf2d.T,'k')
        plot(ratei,np.vstack(pdfi2d).T,'b')
        sca(axL[1])
        plot(rateiconv,pdfconv,'r')
        x = [ratemi]
        y = [1.1*max(pdfconv)]
        xerr = [[-rate_err2],[rate_err1]]
        errorbar(x, y, xerr=xerr, fmt='o')

    out['rate'] = ratemi
    out['rate_err1'] = rate_err1
    out['rate_err2'] = rate_err2
    return out
