import numpy as np
import xarray as xr 

# Some nice, convenient grids
per_bins_dict = {
# sqrt2 
    'xfine': np.round(np.logspace(np.log10(0.1),np.log10(1000),33),4),
# sqrt2 
    'fine': [ 
        1.00, 1.41, 2.00,  2.83,  4.00, 5.66, 8.00,  
        11.3, 16., 22.6, 32.0, 45.3, 64.0, 90.5, 128., 
        181,  256 ],
    'coarse': [ 
        1.00, 2.00,  4.00, 8.00,  
        16., 32.0, 64.0, 128., 256 ],
    'four-per-decade': 10**np.arange(0,3.001,0.25)    
}    

bins = per_bins_dict['xfine']
per_bins_dict['xfine-hj'] = bins[(1 <= bins) & (bins <= 10)]

prad_bins_dict = {
    'xfine': np.round(np.logspace(np.log10(0.5),np.log10(32),49 ),2),
    'fine': np.round(np.logspace(np.log10(0.5),np.log10(32),25 ),2),
    'coarse': [ 
         0.5,  0.71, 1.00, 1.41, 2.0, 2.83, 4.00, 5.66, 8.0, 11.31, 16.0
    ]
}

bins = prad_bins_dict['xfine']
prad_bins_dict['xfine-hj'] = bins[(8 <= bins) & (bins <= 32)]

smet_bins_dict = {
    'fine': np.arange(-0.8,0.6001,0.2),
    'four': np.arange(-0.75,0.451,0.3)
}

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
