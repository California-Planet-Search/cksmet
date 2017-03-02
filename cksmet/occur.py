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
        df (pandas.DataFrame): sample of planets. must contain the following 
            rows:
            - iso_sma (semi-major axis)
            - iso_srad (stellar radius)
            - iso_srad (stellar radius)
            - must also contain quantities associated with the bins. e.g. if 
              we're binning in `period`, must contian `period`
    """

    def __init__(self, plnt):
        self.plnt = plnt

    def add_prob_transit(self, b_transit):
        """Compute transit probability
        
        Args:
            b_transit (float): if we require b < 0.7 then the transit 
                probability is a / Rstar / b_transit
        """

        plnt = self.plnt
        self.b_transit = b_transit
        a = np.array(plnt['iso_sma']) * u.AU
        srad = np.array(plnt['iso_srad']) * c.R_sun
        plnt['prob_transit'] = (srad * b_transit / a).cgs.value
        plnt['inv_prob_transit'] = plnt['prob_transit']**-1
        self.plnt = plnt

    def add_bins(self, bins_dict, spacing_dict):
        """Lay down grid for computing occurrence
        
        Args:
            bins_dict (dict): Dictionary of lists defining bins. e.g.:
                {'period': [5,10,20], 'prad': [1,2,4], 'smet': [-0.4,0.0,0.4]}
            spacing_dict (dict): Specify linear or log spacing e.g.:
                {'period': 'log', 'prad': 'log', 'smet': 'linear'}
        
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

    def count_planets_in_bins(self, sbins):
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
        g = self.plnt.groupby(labels)
        out = g[['id_starname']].count()
        out['plnt_transit_sum'] = g[['inv_prob_transit']].count()
        out['plnt_total_sum'] = g[['inv_prob_transit']].sum()
        out = out.to_xarray()
        for sbin in sbins :
            out.coords[sbin] = self.binsc[sbin]
        return out 

class Stars(object):
    """
    Class for computing the number of stars in my sample around which
    a planet could be detected.
    """

    minimum_mes = 15 # Remove stars that don't achieve mes of 15.
    
    def __init__(self, stars):
        """
        Initialize stellar sample object

        Args:
            stars: DataFrame with the following columns.
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

    def count_stars_in_bins(s):
        """
        Return the number of stars ameable to the detection of a planet of
        a given size, radius, and metallicity
        """
        pass
        
    def f_detectable(self, period, prad):
        """Count number of stars where a planet of a particular size and
        radius could be detected.
        """        
        pass
        
    def mes(self, period, prad):
        """
        Calculate Multiple Event Statistic

        For a planet of a given period and size, calculate the
        expected MES for every star in the sample.

        Args:
            period: orbital period of planet days
            prad: size of planet in Earth-radii 

        Returns:
            array: MES

        """
        depth = self._depth(prad)
        tdur = self._tdur(period)
        cdpp = self._cdpp(tdur) * 1e-6
        num_transits = self._num_transits(period)
        _mes = depth / cdpp * np.sqrt(num_transits)
        return _mes

    def _tdur(self, period):
        """
        Compute duration for a putative transit of a given orbital period
        
        Args:
            period (float): orbital period

        Returns:
            pandas.Series: transit duration for each star in the sample.
        """
        period_yrs = period / 365.25
        tdur = ( 
            TDUR_EARTH_SUN_HRS * self.stars['srad'] * 
            (period_yrs / self.stars['smass'] )**(1.0/3.0)
        )
        return tdur

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

    def _num_transits(self, period):
        """
        Compute number of transits for a planet having a given orbital
        period, factoring in duty cycle.

        Args:
            period (float): orbital period (days)
          
        Returns:
            pd.Series: the number of transits for each star in the sample
        """


        return self.stars.tobs / period

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

