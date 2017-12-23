import cksmet.io
import cksmet.surface
from scipy import ndimage as nd
import numpy as np
from pinky import Pinky
import pandas as pd
import numpy.random
nstars = int(1e5)

def load_pinky():
    df = cksmet.io.load_table('occur-surface',cache=1)
    # convert into an x-array for plotting
    ds = df.groupby(['per1','prad1']).first().to_xarray()
    rate = ds.rate
    rate = rate.fillna(1e-6)

    # Smooth out the contours a bit
    rate = nd.gaussian_filter(rate,(4,2))

    # Only generate planets in the limits of the contour plot
    _m = (np.log10(1) - np.log10(0.5)) / (np.log10(100)-np.log10(1))
    _b = np.log10(0.5)

    b = (
        (1 < ds.perc) & (ds.perc < 350) & 
        (0.5 < ds.pradc) & (ds.pradc < 32) & 
        (np.log10(ds.pradc) > (_m * np.log10(ds.perc) + _b))
    )

    print _m, _b
    rate[~b] = 0.0
    eps = 1e-10
    X, Y = np.log10(ds.perc), np.log10(ds.pradc)

    extent=[float(ds.per1.min()),float(ds.per2.max()),float(ds.prad1.min()),float(ds.prad2.max())]
    extent=np.array(extent)
    extent=np.log10(extent)
    p = Pinky(P=rate.T,extent=extent)
    fac = cksmet.surface.ssamp_per * cksmet.surface.ssamp_prad
    nplanets = int(nstars * rate.sum() / fac)
    return p, nplanets
    
def load_population(p, nplanets):
    np.random.seed(0)
    samples = p.sample(n_samples=nplanets,r=20)
    samples = 10**samples
    samples = pd.DataFrame(samples,columns=['per','prad'])
    
    hotjup = samples.query('1 < per < 10 and 8 < prad < 32')
    print "{} stars ".format(nstars)
    print "{} planets ".format(nplanets)
    print "{:f} hot-Jupiters per 100 stars (should be ~0.57)".format(1.0*len(hotjup)/nstars*100.0)
    return samples

