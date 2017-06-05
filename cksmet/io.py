import os
from collections import Iterable

import pandas as pd
import numpy as np
from scipy.io import idl
from astropy.io import fits 
from astropy import constants as c 
from astropy import units as u
from astropy.table import Table

from cpsutils.pdplus import LittleEndian
import cksphys.io
import cksmet.cuts
import cksspec.io

DATADIR = os.path.join(os.path.dirname(__file__),'../data/')
FLUX_EARTH = (c.L_sun / (4.0 * np.pi * c.au**2)).cgs.value

def load_table(table, cache=0, cachefn='load_table_cache.hdf', verbose=False):
    """Load tables used in cksphys

    Args:
        table (str): name of table. must be one of


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
        
    elif table=='lamost-dr2-cuts':
        # Apply similar set of cuts to the lamost sample.
        df = load_table('lamost-dr2',cache=1)
        cuttypes = cksmet.cuts.lamo_cuttypes
        df = cksmet.cuts.add_cuts(df, cuttypes, 'lamo')

    elif table=='field-cuts':
        # Apply similar set of cuts to the lamost sample.
        df = load_table('huber14+cdpp',cache=1)
        cuttypes = cksmet.cuts.plnt_cuttypes
        df = cksmet.cuts.add_cuts(df, cuttypes, 'field')

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

def split_err(sval):
    """Parse the double-sided errors from hadden""" 
    try:
        istart = sval.index('$')
        istop = sval.index('^')
        val = sval[istart+1:istop]
        sval = sval[istop+1:]

        istart = sval.index('{')
        istop = sval.index('}')
        err1 = sval[istart+1:istop]
        sval = sval[istop+1:]

        istart = sval.index('{')
        istop = sval.index('}')
        err2 = sval[istart+1:istop]
    except ValueError:
        val, err1, err2 = None, None, None
    return val, err1, err2


def table_bin(df, bins):
    g = df.groupby(pd.cut(df.iso_prad,bins=bins))
    g = g['cks_smet']
    dfbin = pd.DataFrame(index=g.first().index)
    dfbin['bin0'] = bins[0:-1]
    dfbin['bin1'] = bins[1:]
    dfbin['binc'] = dfbin.eval('sqrt(bin0*bin1)')
    dfbin['count'] = g.count()
    dfbin['fe_mean'] = g.mean()
    dfbin['fe_std'] = g.std()
    dfbin['fe_mean_err'] = dfbin['fe_std']/ np.sqrt(dfbin['count'])
    dfbin['fe_p05'] = g.apply(lambda x : np.percentile(x, 05))
    dfbin['fe_p95'] = g.apply(lambda x : np.percentile(x, 95))
    dfbin['fe_p01'] = g.apply(lambda x : np.percentile(x, 01))
    dfbin['fe_p99'] = g.apply(lambda x : np.percentile(x, 99))
    dfbin['fe_p50'] = g.apply(lambda x : np.percentile(x, 50))
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



def load_pairs(generate_csv=False):
    pairsfn = 'data/pairs.csv'
    if os.path.exists(pairsfn):
        if not generate_csv:
            print "reading {}".format(pairsfn)
            pairs = pd.read_csv(pairsfn)
            return pairs
    cks = load_table('cks')
    cks.index = cks.kepid
    cks = cks.query('nplanets > 1')
    cks = cks.sort_values(by=['kepid','koi_period'])
    pairs = []
    
    cols_each = 'kepoi_name koi_period iso_prad smax'.split()
    cols_both = 'cks_smet nplanets mass'.split()
    
    irow = 0 
    for kepid in cks.kepid.drop_duplicates():
        nplanets = cks.ix[kepid,'nplanets'].iloc[0]
        i = 0
        while i < nplanets-1:
            row = dict(cks.ix[kepid].iloc[i][cols_both])
            inner = dict(cks.ix[kepid].iloc[i])
            inner['i_planet'] = i+1
            outer = dict(cks.ix[kepid].iloc[i+1])
            outer['i_planet'] = i+2

            for c in cols_each:
                row['inner_'+c] = inner[c]
                row['outer_'+c] = outer[c]

            pairs.append(row)
            i+=1

        if (irow % 100)==0:
            print irow

        irow+=1

    pairs = pd.DataFrame(pairs)
    pairs.to_csv(pairsfn)
    return pairs

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

def load_extended_kic():
    """
    Load KIC parameters with noise information
    """
    #qtbase = np.array(lc_per_quarter.values()) * long_cadence_day


    return df 






