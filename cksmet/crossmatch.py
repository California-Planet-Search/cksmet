"""
Merge LAMOST catalog and K2 catalog
"""

import pandas as pd
import imp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pylab import *
import os
from scipy.spatial import cKDTree

#from terra.utils import plotplus

deg2rad = 2*pi/360
def add_xyz(df):
    # theta = 0 in spherical coords is along z axis
    theta = 90 - df.dec.copy()
    phi = df.ra.copy()

    theta *= deg2rad# [deg] -> radians
    phi *= deg2rad # [deg] -> radians
    df['x'] = sin(theta)*cos(phi)
    df['y'] = sin(theta)*sin(phi)
    df['z'] = cos(theta)
    return df

def match_coords(t0,c0, dist=1.0):
    """
    Match Coordinates

    Given the TESS catalog, find the nearest neighbor on unit sphere
    in the c catalog.  Join the two catalogs together on that nearest
    neighbor. If points difer by more than 1 arcsec drop from
    comparison.

    Args:
        t0 (pandas.DataFrame): Catalog 1. Must contain ra and dec columns
        c0 (pandas.DataFrame): Catalog 2. Must contain ra and dec columns
        dist (float): Distance in arcsc
    """
    for cat in [t0,c0]:
        for col in ['ra','dec']:
            assert list(cat.columns).count(col)==1,"must contain {}".format(col)

    c = c0.copy()
    t = t0.copy()
    c.index=range(len(c))
    t.index=range(len(t)) # i returned from cKDT tree first element in array has index=0
    
    c = add_xyz(c)
    t = add_xyz(t)
    tc = t.copy()
    # Cut out TESS targets that are over 20 degrees center of campagin
    dmed = t['x y z'.split()] - c.median()['x y z'.split()]
    distmed = np.sqrt(dmed['x']**2 + dmed['y']**2 + dmed['z']**2)
    t = t[distmed < 20 * deg2rad]
    
    tree = cKDTree(array(c['x y z'.split()]))
    d,i = tree.query(t['x y z'.split()],distance_upper_bound=dist/206265.)

    tc['cid'] = i # index of EPIC catalog
    tc['cd'] = d*206265. # distance to nearest EPIC target star (arcsec)
    tc = tc[tc.cd < dist]
    tc = pd.merge(tc,c,left_on='cid',right_index=True,suffixes=['_t','_c'])

    pd.set_option('precision',4)


    keys = 'ra dec'.split()
    f = lambda k : [k+'_t',k+'_c'] 
    keys = reduce(lambda x,y : x+y,(map(f,keys)))
    keys += ['cd']

    print "Spot check that coordinates were correctly matched"
    s = tc[keys].head().T.to_string()
    s += " <- displacement between the two catalogs (arcsec)"
    print s
    return tc

def read_tess():
    path = os.path.join(TESS_CATALOGS,'tess_dwarf_star_catalog.csv.gz')
    namemap = {'UID':'uid','RA':'ra','DEC':'dec','I mag':'imag','I mag':'imag',
               'Teff':'teff','V mag':'vmag','J mag':'jmag'}
    cat = pd.read_csv(path,compression='gzip')
    cat = cat.rename(columns=namemap)
    cat = cat[namemap.values()]
    return cat

def sanity_plots(tc):
    fig,axL = subplots(ncols=2,figsize=(8,6))
    sca(axL[0])
    hist(array(tc.cd),bins=100)
    xlabel('Dist from TESS to Nearest EPIC Star (arcsec)')

    sca(axL[1])    
    plot(tc.vmag,tc.kepmag,'.r',alpha=0.2)
    xlabel('TESS Vmag')
    ylabel('EPIC Vmag')

def onSi(tc0_path,flag_path):
    """
    Merge the combined TESS/EPIC lists with the output of the K2FOV plot
    """

    dfonSi = pd.read_csv(tc0_path,index_col=0)
    dfonSi.index=range(len(dfonSi)) # Make sure indecies span [0:len(tc)-1]

    # Read in flag info
    flag = pd.read_csv(flag_path,header=None,names='ra dec kepmag flag'.split())
    flag = flag[['kepmag','flag']]

    dfonSi = pd.merge(dfonSi,flag,left_index=True,right_index=True)
    assert np.allclose(dfonSi.kepmag_x,dfonSi.kepmag_y,atol=0.1),\
        "Indexing error tc kepmag not equal to flag kepmag"
    
    dfonSi = dfonSi.drop('kepmag_y',axis=1)
    dfonSi = dfonSi.rename(columns=dict(kepmag_x='kepmag'))
    return dfonSi

def merge_tess_epic(tess,epic,k2_camp,nlimit=40000):
    """
    Create a CSV file that we can feed into Tom Barclay's runK2onSilicon.py
    """
    figure()
    plot(tess.ra,tess.dec,',')

    figure()
    plot_sample(tess,epic)

    tess_epic = match_coords(tess,epic)

    tess_epic = tess_epic.sort('kepmag')
    print "selecting brightest %i" % nlimit
    tess_epic = tess_epic.iloc[:nlimit]
    kepmag = tess_epic.kepmag
    print "kepmag = %f to %f" % ( min(kepmag), max(kepmag) )

    sanity_plots(tess_epic)

    for suffix,df in zip(['','_debug'],[tess_epic,tess_epic.iloc[::100]]):
        os.system('mkdir -p %s' % k2_camp)
        name = '%s/TESS-Dwarf+EPIC_%s%s.csv' % (k2_camp,k2_camp,suffix)
        df['ra_c dec_c kepmag'.split()].to_csv(name,header=None,index=False)
        print "Just made %s" % name

        name = '%s/TESS-Dwarf+EPIC_%s%s_full.csv' % (k2_camp,k2_camp,suffix)
        df.to_csv(name)
        print "Just made %s" % name


    return tess_epic

def k2_to_excel(dfK2,outfile):
    dfK2.to_excel(outfile.replace('.xls','_full.xls'),'Sheet1')
    dfK2 = dfK2[dfK2.flag.isin([1,2])]
    dfK2 = dfK2['epic kepmag flag'.split()]
    dfK2 = dfK2.sort(['flag','kepmag'],ascending=[False,True])
    dfK2['cad'] = 30
    dfK2['comments'] = dfK2.flag.apply(lambda x : "K2fov_flag=%i" % x)
    dfK2.to_excel(outfile,'Sheet1')

def proposal_plots(k2_camps):
    fig,axL = subplots(
        ncols=len(k2_camps),figsize=(8,3.5)
        )

    for i,k2_camp in enumerate(k2_camps):
        ax = axL[i]
        sca(ax)

        tess_epic_onSi = onSi(
            '%s/TESS-Dwarf+EPIC_%s_full.csv' %(k2_camp,k2_camp),
            '%s/targets_siliconFlag.csv' %(k2_camp) 
             )
        erik = pd.read_excel('K2-%s-Petigura.xls' % k2_camp,'Sheet1')
        tess_epic_onSi.index =tess_epic_onSi.epic
        tess_epic_onSi = tess_epic_onSi.ix[erik.epic]
        tess_epic_onSi = tess_epic_onSi[tess_epic_onSi.flag==2]
        plot(tess_epic_onSi.teff,tess_epic_onSi.kepmag,'.',ms=2,alpha=0.5)    
        setp(ax.xaxis.get_ticklabels(),rotation=40)
        nOnSi = len(tess_epic_onSi)
        plotplus.AddAnchored("%s \nTargets: %i" % (k2_camp,nOnSi),loc=1,
                    prop={'size':'small'} )
        gcf().set_tight_layout(True)
        ylim(4,14)
        plotplus.flip('both')

    axL[0].set_xlabel('Effective Temerature')
    axL[0].set_ylabel('Kepler Magnitude')
    gcf().savefig('ProposalPlots.png',dpi=160)



    
def compute_nearest_neighbors(cat1,cat2,sep):
    """
    cat1 : target catalog
    cat2 : background star catalog
    sep : arcsec

    Returns
    -------
    Rows are every star in cat2 within sep of every star in cat1
    """
    cat1 = cat1.copy()
    cat1.index = range(len(cat1))
    dfnn = []
    kdt = None

    n = 2
    while 1:
        idx1, idx2, ds, kdt = spherematch(
            cat1.ra, cat1.dec, cat2.ra, cat2.dec, tol=None, nnearest=n, kdt=kdt
        )

        temp = pd.DataFrame(cat2['epic ra dec'.split()].ix[idx2])
        temp['ds'] = ds * 60 * 60
        temp['n'] = n
        temp.index = range(len(temp))
        temp = temp[temp.ds < sep]
        
        nstar = len(temp)
        print "%i stars have %i-th nn inside %f arcsec" % (nstar,n,sep)
        if nstar==0:
            break

        test =  pd.merge(
            cat1['epic kepmag ra dec'.split()], temp, left_index=True, 
            right_index=True,suffixes=['','_neighbor']
        )
        dfnn += [test]
        n+=1 
    dfnn = pd.concat(dfnn)
    return dfnn

"""
Match two sets of on-sky coordinates to each other.
I.e., find nearest neighbor of one that's in the other.

Similar in purpose to IDL's spherematch, but totally different implementation.

Requires numpy and scipy.
"""
import numpy as np
try:
    from scipy.spatial import cKDTree as KDT
except ImportError:
    from scipy.spatial import KDTree as KDT
 
 
 
 
def spherematch(ra1, dec1, ra2, dec2, tol=None, nnearest=1,kdt=None):
    """
    Finds matches in one catalog to another.

    Parameters
    ra1 : array-like
        Right Ascension in degrees of the first catalog
    dec1 : array-like
        Declination in degrees of the first catalog (shape of array must match `ra1`)
    ra2 : array-like
        Right Ascension in degrees of the second catalog
    dec2 : array-like
        Declination in degrees of the second catalog (shape of array must match `ra2`)
    tol : float or None, optional
        How close (in degrees) a match has to be to count as a match.  If None,
        all nearest neighbors for the first catalog will be returned.
    nnearest : int, optional
        The nth neighbor to find.  E.g., 1 for the nearest nearby, 2 for the
        second nearest neighbor, etc.  Particularly useful if you want to get
        the nearest *non-self* neighbor of a catalog.  To do this, use:
        ``spherematch(ra, dec, ra, dec, nnearest=2)``

    Returns
    -------
    idx1 : int array
        Indecies into the first catalog of the matches. Will never be
        larger than `ra1`/`dec1`.
    idx2 : int array
        Indecies into the second catalog of the matches. Will never be
        larger than `ra1`/`dec1`.
    ds : float array
        Distance (in degrees) between the matches



    """
 
    ra1 = np.array(ra1, copy=False)
    dec1 = np.array(dec1, copy=False)
    ra2 = np.array(ra2, copy=False)
    dec2 = np.array(dec2, copy=False)
 
    if ra1.shape != dec1.shape:
        raise ValueError('ra1 and dec1 do not match!')
    if ra2.shape != dec2.shape:
        raise ValueError('ra2 and dec2 do not match!')
 
    x1, y1, z1 = _spherical_to_cartesian(ra1.ravel(), dec1.ravel())
 
    # this is equivalent to, but faster than just doing np.array([x1, y1, z1])
    coords1 = np.empty((x1.size, 3))
    coords1[:, 0] = x1
    coords1[:, 1] = y1
    coords1[:, 2] = z1
 
    x2, y2, z2 = _spherical_to_cartesian(ra2.ravel(), dec2.ravel())
 
    # this is equivalent to, but faster than just doing np.array([x1, y1, z1])
    coords2 = np.empty((x2.size, 3))
    coords2[:, 0] = x2
    coords2[:, 1] = y2
    coords2[:, 2] = z2

    if type(kdt)==type(None):
        kdt = KDT(coords2)
    
    if nnearest == 1:
        idxs2 = kdt.query(coords1)[1]
    elif nnearest > 1:
        idxs2 = kdt.query(coords1, nnearest)[1][:, -1]
    else:
        raise ValueError('invalid nnearest ' + str(nnearest))
 
    ds = _great_circle_distance(ra1, dec1, ra2[idxs2], dec2[idxs2])
 
    idxs1 = np.arange(ra1.size)
 
    if tol is not None:
        msk = ds < tol
        idxs1 = idxs1[msk]
        idxs2 = idxs2[msk]
        ds = ds[msk]
 
    return idxs1, idxs2, ds, kdt
 
 
def _spherical_to_cartesian(ra, dec):
    """
    (Private internal function)
    Inputs in degrees.  Outputs x,y,z
    """
    rar = np.radians(ra)
    decr = np.radians(dec)
 
    x = np.cos(rar) * np.cos(decr)
    y = np.sin(rar) * np.cos(decr)
    z = np.sin(decr)
 
    return x, y, z
 
 
def _great_circle_distance(ra1, dec1, ra2, dec2):
    """
    (Private internal function)
    Returns great circle distance.  Inputs in degrees.

    Uses vicenty distance formula - a bit slower than others, but
    numerically stable.
    """
    from numpy import radians, degrees, sin, cos, arctan2, hypot
 
    # terminology from the Vicenty formula - lambda and phi and
    # "standpoint" and "forepoint"
    lambs = radians(ra1)
    phis = radians(dec1)
    lambf = radians(ra2)
    phif = radians(dec2)
 
    dlamb = lambf - lambs
 
    numera = cos(phif) * sin(dlamb)
    numerb = cos(phis) * sin(phif) - sin(phis) * cos(phif) * cos(dlamb)
    numer = hypot(numera, numerb)
    denom = sin(phis) * sin(phif) + cos(phis) * cos(phif) * cos(dlamb)
    return degrees(arctan2(numer, denom))
