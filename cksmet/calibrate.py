"""Derive calibration for SpecMatch parameters. Compare against
touchstone stars to assess the integrity of results.
"""
import os 

import pandas as pd
import cksmet._calibrate
import cksmet.io
import cksmet.utils
from cksmet.io import DATADIR
from matplotlib.pylab import *

CAL_SAMPLE = 'cks' # stellar sample used for calibration
LAMO_SAMPLE = 'lamost-dr2'  # lamost table used for calibration

def add_booleans(df0):
    df = df0.copy()
    print "sm shares {} stars with sx".format(df.teff_sx_uncal.count())
    b = (
        (abs(df.teff_sm - df.teff_sx_uncal) < 300) &
        (abs(df.logg_sm - df.logg_sx_uncal) < 0.35) &
        (abs(df.fe_sm - df.fe_sx_uncal) < 0.3) 
    )
    df['bagree'] = b.astype(float)
    idx = df[df.teff_sx_uncal.isnull()].index
    df.ix[idx,'bagree'] = 2  
    df['bagree'] = df.bagree.astype(int)
    
    df['bteff'] = df.teff_sm.between(4700,6500)
    df['bvsini'] = df.vsini_sm.between(0,20)
    df['bcal'] = (df.bagree==1.0) & df.bteff & df.bvsini 
    return df

def _calibrate_lamo(df, calfn, verbose, fitkw):
    # Perform calibration on teff
    node_points = [
        dict(teff=7000),
        dict(teff=4000)
    ]

    node_points = pd.DataFrame(node_points)
    calteff = cksmet._calibrate.Calibrator('teff')
    calteff.fit(df, node_points, **fitkw)

    # Calibrate logg
    node_points = [
        dict(logg=1.0),
        dict(logg=5.0),
    ]

    node_points = pd.DataFrame(node_points)
    callogg = cksmet._calibrate.Calibrator('logg')
    callogg.fit(df, node_points, **fitkw)

    # Perform calibration on [Fe/H]
    node_points = [
        dict(fe=-1.0),
        dict(fe=1.0),
    ]
    node_points = pd.DataFrame(node_points)
    calfe = cksmet._calibrate.Calibrator('fe')
    calfe.fit(df, node_points, **fitkw)

    # Save away validation parameters
    calteff.to_fits(calfn)
    callogg.to_fits(calfn)
    calfe.to_fits(calfn)

    '''
    # Print calibration planes
    x = dict(teff=5500)
    dx = dict(teff=100)
    calteff.print_hyperplane(x,dx,fmt='.2f')

    x = dict(logg=3.5)
    dx = dict(logg=0.1)
    callogg.print_hyperplane(x,dx,fmt='.4f')

    x = dict(fe=0.0)
    dx = dict(fe=0.1)
    calfe.print_hyperplane(x,dx,fmt='.4f')
    '''

def calibrate_lamo_bootstrap():
    """
    Calibrates LAMOST parameters to CKS scale 
    """
    caltable = 'lamost-cks-calibration-sample'
    calfn_pre  = 'bootstrap/L2_cal_lamo-to-cks'
    df = cksmet.io.load_table(caltable)
    verbose = False

    n = 1000
    np.random.seed(0)
    fitkw=dict(suffixes=['_new','_lib'], mode="L2")
    for isamp in range(n):
        if isamp%10==0:
            print isamp
        df_samp = df.sample(frac=1,replace=True) 
        calfn = "{}_{:04d}.fits".format(calfn_pre,isamp)
        _calibrate_lamo(df_samp, calfn, verbose,fitkw)
    caltable = 'lamost-cks-calibration-sample-noclip'
    calfn_pre  = 'bootstrap/L1_noclip_cal_lamo-to-cks'
    df = cksmet.io.load_table(caltable)
    fitkw=dict(suffixes=['_new','_lib'], mode="L1")
    for isamp in range(n):
        if isamp%10==0:
            print isamp
        df_samp = df.sample(frac=1,replace=True) 
        calfn = "{}_{:04d}.fits".format(calfn_pre,isamp)
        _calibrate_lamo(df_samp, calfn, verbose, fitkw)

def calibrate_lamo():
    """
    Calibrates LAMOST parameters to CKS scale 
    """
    caltable = 'lamost-cks-calibration-sample'
    calfn  = 'cal_lamo-to-cks.fits'

    df = cksmet.io.load_table(caltable)

    verbose = True
    _calibrate_lamo(df, calfn, verbose)
    if verbose:
        print "saving cal file to {}".format(calfn)

    print "Calibrate all of the LAMOST parameters"
    lamo = cksmet.io.load_table(LAMO_SAMPLE,cache=1)
    lamo = cksmet.io.sub_prefix(lamo,'lamo_')
    lamo = lamo['id_kic kic_kepmag steff slogg smet'.split()]
    namemap = {'steff':'teff','slogg':'logg','smet':'fe'}
    lamo = lamo.rename(columns=namemap)
    lamocal = cksmet._calibrate.calibrate(lamo, calfn, mode='uncal')
    lamocal = cksmet.utils.replace_columns(lamocal,'teff','lamo_steff')
    lamocal = cksmet.utils.replace_columns(lamocal,'logg','lamo_slogg')
    lamocal = cksmet.utils.replace_columns(lamocal,'fe','lamo_smet')
    lamocal = lamocal.drop(['delta'],axis=1)

    lamofn = os.path.join(DATADIR,'lamost-dr2-cal.hdf')
    ncal = len(lamocal)
    print "saving {} calibrated LAMO parameters into {}".format(ncal,lamofn)
    lamocal.to_hdf(lamofn,'lamost-dr2-cal')
