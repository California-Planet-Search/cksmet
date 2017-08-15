import os
import cPickle as pickle

import pandas as pd
import numpy as np
from scipy.io import idl
from astropy.io import fits 
from astropy import constants as c 
from astropy import units as u

import cksmet.cuts
import cksphys.io
import cksspec.io
import cksmet.grid
import cksmet.analysis

DATADIR = os.path.join(os.path.dirname(__file__),'../data/')
FLUX_EARTH = (c.L_sun / (4.0 * np.pi * c.au**2)).cgs.value
COMPLETENESS_FILE = os.path.join(DATADIR, 'comp.pkl')

def load_table(table, cache=0, cachefn='load_table_cache.hdf', verbose=False):
    """Load tables used in cksphys

    Args:
        table (str): name of table. must be one of
            - nea 


        cache (Optional[int]): whether or not to use the cache
            - 0: don't use the cache recreate all files
            - 1: read from cache
            - 2: write tables to cache

    Returns:
        pandas.DataFrame: table

    """
    if cache==1:
        try:
            df = pd.read_hdf(cachefn,table)
            print "read table {} from {}".format(table,cachefn)
            return df
        except IOError:
            print "Could not find cache file: %s" % cachefn
            print "Building cache..."
            cache=2
        except KeyError:
            print "Cache not built for table: %s" % table
            print "Building cache..."
            cache=2

    if cache==2:
        df = load_table(table, cache=False)
        print "writing table {} to cache".format(table)
        df.to_hdf(cachefn,table)
        return df

    if table=='nea':
        df = cksphys.io.load_table('nea')

    elif table=='cksbin-fe':
        bins = [0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 5.7, 8.0, 11.3, 16]
        cks = load_table('cks')
        #cks = cks.query(' -0.3 < feh_cks < 0.5')
        cksbin = table_bin(cks, bins)
        df = cksbin

    elif table=='cksbin-nplanets':
        cks = load_table('cks')
        g = cks.groupby('nplanets')
        g = g['feh_cks']
        dfbin = pd.DataFrame(index=g.first().index)
        print dfbin
        dfbin['multiplicity'] = dfbin.index
        dfbin['nplanets'] = g.count() 
        dfbin['nstars'] = dfbin['nplanets']/ dfbin.multiplicity
        dfbin['fe_mean'] = g.mean()
        dfbin['fe_std'] = g.std()
        dfbin['fe_mean_err'] = dfbin['fe_std']/ np.sqrt(dfbin['nstars'])
        df = dfbin

    elif table=='lamost-dr2':
        df = cksspec.io.load_table('lamost-dr2')

    elif table=='lamost-dr2-cal':
        fn = os.path.join(DATADIR,'lamost-dr2-cal.hdf')
        df = pd.read_hdf(fn, 'lamost-dr2-cal')

    elif table=='cks':
        df = pd.read_csv('../CKS-Physical/data/cks_physical_merged.csv')

    elif table=='lamost-dr2+cks':
        cks = load_table('cks')
        lamost = load_table('lamost-dr2',cache=1)
        df = pd.merge(lamost,cks,on='id_kic')
        return df

    elif table=='cks-cuts':
        df =  cksphys.io.load_table(
            'cks+nea+iso-floor+huber-phot+furlan',cache=1,
            cachefn='../CKS-Physical/load_table_cache.hdf'
        )
        cuttypes = 'none notdwarf faint diluted grazing longper badprad'.split()
        df = cksmet.cuts.add_cuts(df, cksmet.cuts.plnt_cuttypes, 'cks')
        
    elif table=='lamost-dr2-cal-cuts':
        # Apply similar set of cuts to the lamost sample.
        df = load_table('lamost-dr2-cal')
        cuttypes = cksmet.cuts.lamo_cuttypes
        df = cksmet.cuts.add_cuts(df, cuttypes, 'lamo')


    elif table=='lamost-dr2-cal-cuts+cdpp':
        lamo = cksmet.io.load_table('lamost-dr2-cal-cuts',cache=1)
        huber14 = cksmet.io.load_table('huber14+cdpp',cache=1)
        huber14 = huber14['id_kic kepmag cdpp3 cdpp6 cdpp12'.split()]
        lamo = pd.merge(lamo,huber14,on='id_kic')
        df = lamo

    elif table=='field-cuts':
        # Apply similar set of cuts to the KIC sample
        df = load_table('huber14+cdpp',cache=1)
        cuttypes = cksmet.cuts.plnt_cuttypes
        df = cksmet.cuts.add_cuts(df, cuttypes, 'field')

    elif table=='cdpp':
        fn = os.path.join(DATADIR,'kic_q0_q17.dat')
        df = idl.readsav(fn)
        df = pd.DataFrame(df['kic'])

    elif table=='huber14+cdpp':
        df = load_table('cdpp',cache=1)
        df = df.rename(columns={
            'KEPMAG':'kepmag','KICID':'id_kic','CDPP3':'cdpp3','CDPP6':'cdpp6',
            'CDPP12':'cdpp12'}
        )
        df = df['id_kic kepmag cdpp3 cdpp6 cdpp12'.split()]
        cdpp = np.vstack(df.ix[:,'cdpp3'])
        for col in 'cdpp3 cdpp6 cdpp12'.split():
            cdpp = np.vstack(df.ix[:,col])
            cdpp[cdpp==0.0] = np.nan
            cdpp = np.nanmedian(cdpp,axis=1)
            df[col] = cdpp
            df['log'+col] = np.log10(cdpp)

        huber14 = cksspec.io.load_table('huber14')
        stellar17 = cksspec.io.load_table('stellar17')
        df = pd.merge(df,huber14)
        df = pd.merge(df,stellar17)
        
        # add place holder column for quarter zero
        df['st_quarters'] = '0'+df.st_quarters 
        df['tobs'] = 0
        for q in range(18):
            # Was quarter observed
            qobs = df.st_quarters.str.slice(start=q,stop=q+1).astype(int)
            if q==17:
                qobs=0
            df['tobs'] += qobs * long_cadence_day * lc_per_quarter[q]

    else:
        assert False, "table {} not valid table name".format(table)
    return df


def add_prefix(df,prefix,ignore=['id']):
    namemap = {}
    for col in list(df.columns):
        skip=False
        for _ignore in ignore:
            if col.count(_ignore) > 0:
                skip = True
        if not skip:
            namemap[col] = prefix + col 
    df = df.rename(columns=namemap)
    return df

def sub_prefix(df, prefix,ignore=['id']):
    namemap = {}
    for col in list(df.columns):
        skip=False
        for _ignore in ignore:
            if col.count(_ignore) > 0:
                skip = True
        if not skip:
            namemap[col] = col.replace(prefix,'') 
    df = df.rename(columns=namemap)
    return df

def load_object(key, cache=1):
    pklfn = os.path.join(DATADIR,key+'.pkl')
    if cache==1:
        with open(pklfn,'r') as f:
            obj = pickle.load(f)
            return obj

    elif cache==2:
        obj = load_object(key,cache=0)
        with open(pklfn,'w') as f:
            pickle.dump(obj,f)
        
        return obj

    if key.count('occur')==1:
        obj = cksmet.analysis.load_occur(key)
    elif key.count('comp')==1:
        obj = cksmet.analysis.load_completeness()
    elif key.count('fit')==1:
        obj = cksmet.analysis.load_fit(key)
    return obj

def load_occur(key, cache=1):
    pklfn = os.path.join(DATADIR,key+'.pkl')
    if cache==1:
        with open(pklfn,'r') as f:
            fit = pickle.load(f)
            return fit 

    elif cache==2:
        fit = load_occur(key,cache=0)
        with open(pklfn,'w') as f:
            pickle.dump(fit,f)


def table_bin(df, bins, key):
    g = df.groupby(pd.cut(df.iso_prad,bins=bins))
    g = g[key]
    dfbin = pd.DataFrame(index=g.first().index)
    dfbin['bin0'] = bins[0:-1]
    dfbin['bin1'] = bins[1:]
    dfbin['binc'] = dfbin.eval('sqrt(bin0*bin1)')
    dfbin['count'] = g.count()
    dfbin[key + '_mean'] = g.mean()
    dfbin[key + '_std'] = g.std()
    dfbin[key +'_mean_err'] = dfbin[key + '_std']/ np.sqrt(dfbin['count'])

    for p in [1,5,25,50,75,95,99]:
        k = key+'_{:02d}'.format(p)
        dfbin[k] = g.quantile(q=p * 0.01)
    return dfbin

def latex_table(table):
    tab = load_table(table)
    if table=='cksbin-fe':
        for i, row in tab.iterrows():
            line = "{bin0:.1f}--{bin1:.1f} &"
            line+=" {count:.0f} & "
            line+=" {fe_mean:+.3f} & "
            line+=" {fe_mean_err:.3f} &"
            line+=" {fe_p01:+.3f} &"
            line+=" {fe_p50:+.3f} &"
            line+=" {fe_p99:+.3f}"
            line+=r"\\"
            line = line.format(**row)
            print line

# Pulled from Kepler data characteristics
lc_per_quarter = {
    0:476,
    1:1639,
    2:4354,
    3:4370,
    4:4397,
    5:4634,
    6:4398,
    7:4375,
    8:3279,
    9:4768,
    10:4573,
    11:4754,
    12:4044,
    13:4421,
    14:4757,
    15:4780,
    16:4203,
    17:1556,
}
long_cadence_day = 29.7 / 60.0 / 24.0 # long cadence measurement in days



