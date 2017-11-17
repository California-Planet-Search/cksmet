"""Derive calibration for SpecMatch parameters. Compare against
touchstone stars to assess the integrity of results.
"""
import os 

import pandas as pd
import smsyn.calibrate
import cksspec.plotting.catalogs
import cksspec.io
import cksmet.io
import cksspec.utils
from cksmet.io import DATADIR
from matplotlib.pylab import *

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

def calibrate_lamo():
    """
    Calibrates LAMOST parameters to CKS scale 

    - name
    - obs
    - teff_sm_uncal
    - logg_sm_uncal
    - fe_sm_uncal
    - vsini_sm
    - delta
    - teff_sm,
    - logg_sm
    - fe_sm
    """
    df = cksmet.io.load_table('lamost-cks-calibration-sample')
    import pdb;pdb.set_trace()
    calfn  = 'cal_lamo-to-cks.fits'
    print "saving cal file to {}".format(calfn)

    #df['fe_lib'] = df['fe_new']

    # Perform calibration on teff
    node_points = [
        dict(teff=7000),
        dict(teff=4000)
    ]

    node_points = pd.DataFrame(node_points)
    calteff = smsyn.calibrate.Calibrator('teff')
    calteff.fit(df, node_points)

    # Calibrate logg
    node_points = [
        dict(logg=1.0),
        dict(logg=5.0),
    ]

    node_points = pd.DataFrame(node_points)
    callogg = smsyn.calibrate.Calibrator('logg')
    callogg.fit(df, node_points)

    # Perform calibration on [Fe/H]
    node_points = [
        dict(fe=-1.0),
        dict(fe=1.0),
    ]
    node_points = pd.DataFrame(node_points)
    calfe = smsyn.calibrate.Calibrator('fe')
    calfe.fit(df, node_points)

    # Save away validation parameters
    calteff.to_fits(calfn)
    callogg.to_fits(calfn)
    calfe.to_fits(calfn)

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
    
    # lamocal = lamo['id_kic teff logg fe'.split()]
    print "Calibrate all of the LAMOST parameters"

    import pdb;pdb.set_trace()
    lamo = cksmet.io.load_table(LAMO_SAMPLE,cache=1)
    lamo = cksmet.io.sub_prefix(lamo,'lamo_')
    lamo = lamo['id_kic kic_kepmag steff slogg smet'.split()]
    namemap = {'steff':'teff','slogg':'logg','smet':'fe'}
    lamo = lamo.rename(columns=namemap)
    lamocal = smsyn.calibrate.calibrate(lamo, calfn, mode='uncal')
    lamocal = cksspec.utils.replace_columns(lamocal,'teff','lamo_steff')
    lamocal = cksspec.utils.replace_columns(lamocal,'logg','lamo_slogg')
    lamocal = cksspec.utils.replace_columns(lamocal,'fe','lamo_smet')
    lamocal = lamocal.drop(['delta'],axis=1)

    lamofn = os.path.join(DATADIR,'lamost-dr2-cal.hdf')
    ncal = len(lamocal)
    print "saving {} calibrated LAMO parameters into {}".format(ncal,lamofn)
    lamocal.to_hdf(lamofn,'lamost-dr2-cal')


