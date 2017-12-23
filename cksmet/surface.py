"""
module to compute occurrence surface
"""
import numpy as np
import cksmet.analysis 
import pandas as pd

bin_per = 0.25 # period bin width in dex
bin_prad = 0.1 # prad bin width in dex
ssamp_per = 10 # supersample factor
ssamp_prad = 10 # supersample factor
per1 = 0.5
per2 = 1000
prad1 = 0.25
prad2 = 32
eps = 1e-3
smet_bins = [-1,0.5]

def compute_occur_surface():
    """
    Compute occurrence surface

    Used for contour plots and samples
    """

    df = []
    count = 0 
    for i in range(ssamp_per):
        for j in range(ssamp_prad):
            if count < 2:
                verbose=1
            else:
                verbose=0

            shift_logper = bin_per / ssamp_per * i
            shift_logprad = bin_prad / ssamp_prad * j
            logper1 = np.log10(per1) + shift_logper
            logper2 = np.log10(per2) + shift_logper
            logprad1 = np.log10(prad1) + shift_logprad
            logprad2 = np.log10(prad2) + shift_logprad
            per_bins = 10**(np.arange(logper1,logper2+eps,bin_per))
            prad_bins = 10**(np.arange(logprad1,logprad2+eps,bin_prad))
            print per_bins
            print prad_bins

            occ = cksmet.analysis.compute_occur_bins(
                per_bins, prad_bins, smet_bins, verbose=verbose
            )
            occ.df['count'] = count
            df.append(occ.df)
            count+=1
            
    df = pd.concat(df)
    return df

