import numpy as np
from astropy import units as u
from astropy import constants as c
import pandas as pd
import xarray as xr
from scipy.interpolate import RegularGridInterpolator

TDUR_EARTH_SUN_HRS = (
    ((4 * c.R_sun**3 * 1.0*u.yr / np.pi / c.G / (1.0*c.M_sun))**(1.0/3.0)).to(u.hr)).value

DEPTH_EARTH_SUN = ((c.R_earth / c.R_sun)**2).cgs.value

per_bins_dict = {
    'xfine': [
        1.  ,    1.19,    1.41,    1.68,    2.  ,    2.38,    2.83,
        3.36,    4.  ,    4.76,    5.66,    6.73,    8.  ,    9.51,
        11.31,   13.45,   16.  ,   19.03,   22.63,   26.91,   32.  ,
        38.05,   45.25,   53.82,   64.  ,   76.11,   90.51,  107.63,
        128.  ,  152.22,  181.02,  215.27,  256.  
    ],
    'fine': [ 
        1.00, 1.41, 2.00,  2.83,  4.00, 5.66, 8.00,  
        11.3, 16., 22.6, 32.0, 45.3, 64.0, 90.5, 128., 
        181,  256 ],
    'coarse': [ 
        1.00, 2.00,  4.00, 8.00,  
        16., 32.0, 64.0, 128., 256 ],
}    

prad_bins_dict = {
    'xxfine': [
         0.5  ,   0.522,   0.545,   0.569,   0.595,   0.621,   0.648,
         0.677,   0.707,   0.738,   0.771,   0.805,   0.841,   0.878,
         0.917,   0.958,   1.   ,   1.044,   1.091,   1.139,   1.189,
         1.242,   1.297,   1.354,   1.414,   1.477,   1.542,   1.61 ,
         1.682,   1.756,   1.834,   1.915,   2.   ,   2.089,   2.181,
         2.278,   2.378,   2.484,   2.594,   2.709,   2.828,   2.954,
         3.084,   3.221,   3.364,   3.513,   3.668,   3.83 ,   4.   ,
         4.177,   4.362,   4.555,   4.757,   4.967,   5.187,   5.417,
         5.657,   5.907,   6.169,   6.442,   6.727,   7.025,   7.336,
         7.661,   8.   ,   8.354,   8.724,   9.11 ,   9.514,   9.935,
        10.375,  10.834,  11.314,  11.815,  12.338,  12.884,  13.454,
        14.05 ,  14.672,  15.322,  16.   ],
    'xfine': [ 
        0.5 ,   0.55,   0.59,   0.65,   0.71,   0.77,   0.84,   0.92,
        1.  ,   1.09,   1.19,   1.3 ,   1.41,   1.54,   1.68,   1.83,
        2.  ,   2.18,   2.38,   2.59,   2.83,   3.08,   3.36,   3.67,
        4.  ,   4.36,   4.76,   5.19,   5.66,   6.17,   6.73,   7.34,
        8.  ,   8.72,   9.51,  10.37,  11.31,  12.34,  13.45,  14.67,  16.
    ],
    'fine': [ 
         0.5, 0.59,   0.71, 0.84, 1.00, 1.19, 1.41, 1.68,
         2.0, 2.38,   2.83, 3.36, 4.00, 4.76, 5.66, 6.73,
         8.0, 9.51,  11.31, 13.45, 16.0],
    'coarse': [ 
         0.5,  0.71, 1.00, 1.41, 2.0, 2.83, 4.00, 5.66, 8.0, 11.31, 16.0
    ]
}

smet_bins_dict = {
    'fine': np.arange(-0.8,0.6001,0.2)
}

def _sbins_to_sbins_id(sbins):
    """Convenience function for converting bin names to bins_id names"""
    return ['{}_bin_id'.format(sbin) for sbin in sbins]


class Occurrence():
    """
    Occurrence Object

    This object allows one to compute planet occurrence in bins of
    orbital period, planet radius, and stellar metallicity

    Args:
        df (pandas.DataFrame): Planets. Must contain the following rows:
            - srad: stellar radius (Solar-radii)
            - per: planet orbital period (days)
            - prad: planet radius (Earth-radii) 
            - smax: planet semi-major axis (AU)
            - mes: planet multiple event statistic
            - impact: planet impact parameter
        prob_det_mes_name (str): Probability of detection as function of MES
        imppar_transit (float): Maximum impact parameter allowed.

    """
    plnt_columns = ['srad','per','prad','smax','impact','mes']
    prob_det_mes_names = ['step-12','step-15']

    def __init__(self, plnt, prob_det_mes_name, impact_transit):
        for col in self.plnt_columns:
            assert list(plnt).count(col)==1, \
                "DataFrame must contain {}".format(col)

        assert self.prob_det_mes_names.count(prob_det_mes_name)==1,\
            "prob_det_mes_name must be one of\n" + "\n".join(prob_det_mes_names)

        idxdrop = plnt[plnt.impact > impact_transit].index
        ndrop = len(idxdrop)
        print "dropping {} planets, impact > {}".format(ndrop,impact_transit)

        plnt = plnt.drop(idxdrop)
        if prob_det_mes_name=='step-12':
            idxdrop = plnt[plnt.mes < 12].index
            ndrop = len(idxdrop)
            print "dropping {} planets, mes < 12".format(ndrop,12)
        elif prob_det_mes_name=='step-15':
            idxdrop = plnt[plnt.mes < 15].index
            ndrop = len(idxdrop)
            print "dropping {} planets, mes < 15".format(ndrop,12)
        
        plnt = plnt.drop(idxdrop)

        self.plnt = plnt
        self.prob_det_mes_name = prob_det_mes_name
        self.impact_transit = impact_transit


    def add_bins(self, bins_dict, spacing_dict):
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
            self._label_bins(key)
            
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

    def _label_bins(self, key):
        """
        Assign rows to bins

        Given a columnn `key` and `bins` figure out which bin each row
        belongs to. Also for every row, return the edges and center of the
        bins.

        Args:
            df (DataFrame): Must contain `key` as a column
            key (str): string label of bin
            bins (array): list of bin edges 
            spacing (str): linear/log how to calculate bin center

        """
        sbin_id = '{}_bin_id'.format(key)
        sbin1 = '{}_bin1'.format(key)
        sbin2 = '{}_bin2'.format(key)
        sbinc = '{}_binc'.format(key)
        plnt = self.plnt
        for col in [sbin_id, sbin1, sbin2, sbinc]:
            plnt[col] = None

        bins = self.bins[key]
        for i in range(len(bins)-1):
            bin1 = self.bins1[key][i]
            bin2 = self.bins2[key][i]
            binc = self.binsc[key][i]
            idx = plnt[plnt[key].between(bin1,bin2)].index
            plnt.ix[idx,sbin_id] = i
            plnt.ix[idx,sbin1] = bin1
            plnt.ix[idx,sbin2] = bin2
            plnt.ix[idx,sbinc] = binc
 
        self.plnt = plnt

    def add_prob(self, st):
        """For each planet compute transit and detection probability
        """
        plnt = self.plnt
        a = np.array(plnt['smax']) * u.AU
        srad = np.array(plnt['srad']) * c.R_sun
        plnt['prob_tr'] = (srad * self.impact_transit / a).cgs.value
        plnt['inv_prob_tr'] = plnt['prob_tr']**-1

        i = 0
        for i, row in plnt.iterrows():
            plnt.ix[i,'prob_det'] = st.prob_det(row.per,row.prad)
            if i%100==0:
                print row.per,row.prad,st.prob_det(row.per,row.prad)

        plnt['prob_total'] = plnt['prob_tr'] * plnt['prob_det']
        for suffix in 'prob_tr prob_det prob_total'.split():
            plnt['inv_'+suffix] = plnt[suffix]**-1.0
        self.plnt = plnt

    def compute_occurrence(self, sbins, st):
        """Count planets in bins

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
        labels = [pd.cut(self.plnt[sbin],self.bins[sbin]) for sbin in sbins]
        g = self.plnt.groupby(labels,as_index=True)
        grid = g[['id_starname']].count()
        
        # Transiting and detected planets per bin
        grid['plnt_trdet_sum'] = g[['inv_prob_transit']].sum()

        # Transiting planets per bin (accounting for missed transiting planets)
        grid['plnt_tr_sum'] = g[['inv_prob_det']].sum()

        # Planets per bin (accounting for missed and non-transiting planets)
        grid['plnt_sum'] = g[['inv_prob_total']].sum()
        grid['plnt_occur'] = grid['plnt_sum'] / st.nstars

        grid = grid.to_xarray()
        for sbin in sbins :
            grid.coords[sbin] = self.binsc[sbin]

        for key in self.bins.keys():
            coord = {key: self.binsc[key]}
            grid[key+'c'] = xr.DataArray(self.binsc[key],coord)
            grid[key+'1'] = xr.DataArray(self.bins1[key],coord)
            grid[key+'2'] = xr.DataArray(self.bins2[key],coord)
        self.grid = grid

    def prob_det_grid(self, st):
        """Adds a field to the grid attribute
        """
        grid = self.grid
        grid = grid.to_dataframe()
        grid['prob_det'] = -1 

        i = 0
        print "perc pradc prob_det"
        for idx, row in grid.iterrows():
            grid.ix[idx,'prob_det'] = st.prob_det(row.perc, row.pradc)
            s ="{perc:.3f} {pradc:.3f} {prob_det:.3f}".format(
                **grid.ix[idx] )
            if i%100==0:
                print s
            i+=1

        self.grid = grid.to_xarray()

    def prob_tr_grid(self, st):
        """Adds a field to the grid attribute
        """
        grid = self.grid
        grid = grid.to_dataframe()
        grid['prob_tr'] = -1 

        i = 0
        s ="perc pradc prob_tr"
        print s
        for idx,row in grid.iterrows():
            grid.ix[idx,'prob_tr'] = st.prob_tr(row.perc)
            s ="{perc:.3f} {pradc:.3f} {prob_tr:.3f}".format(
                **grid.ix[idx] )
            if i%100==0 :
                print s
            i+=1

        self.grid = grid.to_xarray()


class Stars(object):
    """
    Class for computing the number of stars in my sample around which
    a planet could be detected.
    """
    
    def __init__(self, stars, occur):
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
        self.occur = occur
        self.nstars = len(stars)

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
        if self.occur.prob_det_mes_name=='step-15':
            out[mes > 15] = 1.0
        elif self.occur.prob_det_mes_name=='step-12':
            out[mes > 12] = 1.0

        return out

    def prob_det(self, per, prad):
        """Probability that a planet would be detectable

        Probability that transiting planet with orbital period `per`
        and planet radius `prad` would be detectable around a
        randomly-drawn star from the stellar sample.

        Args:
            per (float): orbital period
            prad (float): planet size

        Returns:
            float: fraction of stars in sample where we could have detected 
                planet

        """        
        mes = self.mes_scaled(per, prad)
        _prob_det = 1.0*self.prob_det_mes(mes).sum() / self.nstars 
        return _prob_det

    def prob_tr(self, per):
        a = self._smax(per)
        srad = (np.array(self.stars['srad']) * u.R_sun).to(u.AU).value
        _prob_tr = srad * self.occur.impact_transit / a
        return _prob_tr.mean()

    def compute_mes_factor(self,plnt_mes):
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
        plnt_mes['mes_formula'] = np.nan

        i = 0
        for id_koicand, row in plnt_mes.iterrows():
            try:
                mes_formula = self.mes(row.per,row.prad)
                mes_formula = mes_formula.ix[row.id_kic]
                mes_pipeline = row.mes
                plnt_mes.ix[id_koicand,'mes_formula'] = mes_formula
                plnt_mes.ix[id_koicand,'mes_pipeline'] = mes_pipeline 
                if i< 20:
                    s =  "{} {:.2f} {:.2f}".format(
                        id_koicand, mes_pipeline, mes_formula
                    )
                    print s
            except KeyError:
                pass
        
            i+=1

        self.plnt_mes = plnt_mes
        logmes_pipeline = np.log10(plnt_mes.mes_pipeline)
        logmes_formula = np.log10(plnt_mes.mes_formula)
        logmes_diff_med = np.nanmedian(logmes_formula - logmes_pipeline)
        logmes_diff_std = np.nanstd(logmes_formula - logmes_pipeline)
        print "med(log(mes_formula/mes_pipeline)) {:.2f} (dex)".format(
            logmes_diff_med
        )
        self.logmes_diff_std = logmes_diff_med
        self.logmes_diff_med = logmes_diff_std 
        self.mesfac = 10**(-1.0 * logmes_diff_med)

        i = 0
        for id_koicand, row in plnt_mes.iterrows():
            try:
                mes_formula = self.mes_scaled(row.per,row.prad)
                mes_formula = mes_formula.ix[row.id_kic]
                mes_pipeline = row.mes
                plnt_mes.ix[id_koicand,'mes_formula_scaled'] = mes_formula
                plnt_mes.ix[id_koicand,'mes_pipeline'] = mes_pipeline 
                if i< 20:
                    s =  "{} {:.2f} {:.2f}".format(
                        id_koicand, mes_pipeline, mes_formula
                    )
                    print s
            except KeyError:
                pass
        
            i+=1
        
                
