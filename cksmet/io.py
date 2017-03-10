import pandas as pd
from astropy import constants as c 
from astropy import units as u
import numpy as np
from collections import Iterable
FLUX_EARTH = (c.L_sun / (4.0 * np.pi * c.au**2)).cgs.value
import os
import cksphys.io
from astropy.table import Table
from astropy.io import fits 
from cpsutils.pdplus import LittleEndian
import cksmet.cuts
from scipy.io import idl

def load_table(table, verbose=False, cache=False):
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

    hdffn = 'load_table_cache.hdf'
    if cache==1:
        df = pd.read_hdf(hdffn,table)
        print "read table {} from {}".format(table,hdffn)
        return df
    if cache==2:
        df = load_table(table, cache=False)
        print "writing table {} to cache".format(table)
        df.to_hdf(hdffn,table)
        return df

    '''
    elif table=='nea':
        df = pd.read_csv('data/nea_2017-01-28.csv',comment='#')
    '''
    if table=='nea':
        df = cksphys.io.load_table('nea')

    elif table=='hadden-ttv':
        tablefn = 'data/hadden17/masstable.tex'
        with open(tablefn) as f:
            lines = f.readlines()

        df = []
        for row in lines:
            d = {}
            if row.count('&')!=8:
                continue
            row = row.split('&')
            name = row[0]
            if name.count('*')==0:
                continue

            name = name.split(' ')
            name = name[0]+' '+name[1][0] 
            d['pl_name'] = name
            d['pl_orbper'] = row[1]
            d['pl_rade'],d['pl_radeerr1'],d['pl_radeerr2'] = split_err(row[2])
            d['st_mass'],d['st_masserr1'],d['st_masserr2'] = split_err(row[3])
            d['pl_masse'],d['pl_masseerr1'],d['pl_masseerr2'] = split_err(row[4])
            df.append(d)
        df = pd.DataFrame(df)
        df = df.convert_objects(convert_numeric=True)

    elif table=='hadden-rv':
        tablefn = 'data/hadden17/rvtable.tex'
        with open(tablefn) as f:
            lines = f.readlines()

        df = []
        for row in lines:
            d = {}
            if row.count('&')!=4:
                continue
            row = row.split('&')
            name = row[0]

            name = name.split(' ')
            name = name[0]+' '+name[1][0] 
            d['pl_name'] = name
            d['pl_orbper'] = row[1]
            d['pl_masse'],d['pl_masseerr1'],d['pl_masseerr2'] = split_err(row[2])
            d['pl_rade'],d['pl_radeerr1'],d['pl_radeerr2'] = split_err(row[3])
            df.append(d)
        df = pd.DataFrame(df)
        df = df.convert_objects(convert_numeric=True)

    elif table=='hadden':
        ttv = load_table('hadden-ttv')
        rv = load_table('hadden-rv')
        ttv['pl_massmeth'] = 'ttv'
        rv['pl_massmeth'] = 'rv'
        df = pd.concat([ttv,rv])

    elif table=='hadden+cks':
        # Return hadden TTV planets with Kepler and KOI designations
        df = load_table('hadden')
        cks = cksphys.io.load_table('cksphys-merged',cache=2)
        df1 = pd.merge(
            df, cks, left_on='pl_name', right_on='id_kepler_name',how='left'
        )
        df = pd.concat([df1])

    elif table=='hadden-ttv+cks':
        # Return hadden TTV planets with Kepler and KOI designations
        df = load_table('hadden-ttv')
        cks = cksphys.io.load_table('cksphys-merged',cache=2)
        df1 = pd.merge(
            df, cks, left_on='pl_name', right_on='id_kepler_name',how='left'
        )
        df = pd.concat([df1])

    elif table=='hadden-rv+cks':
        # Return hadden TTV planets with Kepler and KOI designations
        df = load_table('hadden-rv')
        cks = cksphys.io.load_table('cksphys-merged',cache=2)
        df1 = pd.merge(
            df, cks, left_on='pl_name', right_on='id_kepler_name',how='left'
        )
        df = pd.concat([df1])

    elif table=='hadden+cks+nea':
        df = load_table('hadden+cks')
        nea = load_table('nea')
        df.index = df.pl_name
        nea.index = nea.pl_name
        df['st_metfe'] = df['cks_smet']
        df['st_teff'] = df['cks_steff']
        df['st_metfe'] = df.st_metfe.fillna(value=nea.st_metfe)
        df['st_teff'] = df.st_metfe.fillna(value=nea.st_teff)
        df['pl_name st_teff st_metfe'.split()]
        return df

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
        hduL = fits.open('data/lamost/dr2_stellar.fits.gz')
        df = pd.DataFrame(LittleEndian(hduL[1].data))
        df = df.query('tsource=="Kepler"')
        df = df[df.tcomment.str.contains('kplr\d{9}')]
        df['id_kic'] = df.tcomment.str.slice(start=4).astype(int)
        df = df['id_kic teff teff_err logg logg_err feh feh_err'.split()]
        df['steff'] = df['teff']
        df['slogg'] = df['logg']
        df['smet'] = df['feh']
        df['steff_err1'] = df['teff_err']
        df['steff_err2'] = -1.0 * df['teff_err']
        df['slogg_err1'] = df['logg_err']
        df['slogg_err2'] = -1.0 * df['logg_err']
        df['slogg_err1'] = df['logg_err']
        df['slogg_err2'] = -1.0 * df['logg_err']
        df['smet_err1'] = df['feh_err']
        df['smet_err2'] = -1.0 * df['feh_err']
        df = add_prefix(df,'lamo_')

        stellar = cksphys.io.load_table('stellar17',cache=True)
        df = pd.merge(df,stellar)
        return df

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
        df = cksmet.cuts.add_cuts(df,cuttypes)
        
    elif table=='lamost-dr2-cuts':
        # Apply similar set of cuts to the lamost sample.
        df = load_table('lamost-dr2',cache=1)
        cuttypes = cksmet.cuts.lamo_cuttypes
        df = cksmet.cuts.add_cuts(df, cuttypes, 'lamo')

    elif table=='buch14':
        df = pd.read_table(
            'data/buchhave_2014_nature_table-s1.txt',sep='\s*',
            names=[
                'id_koicand','steff','steff_err','slogg','slogg_err','smet',
                'smet_err','smass','smass_err1','smass_err2','srad',
                'srad_err1','srad_err2','prad','prad_err'
            ],
            comment='#'
        )
        df = add_prefix(df,'spc_')
    elif table=='buch14-stars':
        df = load_table('buch14')
        df['id_koi'] = df.id_koicand.str.slice(start=1,stop=6).astype(int)
        df = df.drop(['id_koicand','spc_prad','spc_prad_err'],axis=1)
        df = df.drop_duplicates()
    elif table=='buch14-stars+cks':
        df = load_table('buch14')
        cks = load_table('cks')
        df = pd.merge(df,cks)
    elif table=='kic':
        df = pd.read_hdf('data/kic.hdf','kic')
        namemap = {'KICID':'id_kic','TEFF':'steff','LOGG':'slogg','FEH':'smet'}
        df = df.rename(columns=namemap)[namemap.values()]
        df = add_prefix(df,'kic_')

    elif table=='kic+cks':
        df = load_table('kic')
        cks = load_table('cks')
        df = pd.merge(df, cks)
    elif table=='huber14':
        t = Table.read(
            'data/huber14/table5.dat',
            readme='data/huber14/ReadMe',
            format='ascii.cds'
        )
        namemap = {
            'KIC':'id_kic','Teff1':'steff','e_Teff1':'steff_err1',
            'log.g1':'slogg','e_log.g1':'slogg_err',
            '[Fe/H]1':'smet','e_[Fe/H]1':'smet_err',
            'R':'srad','e_R':'srad_err','M':'smass','e_M':'smass_err'
        }
        df = t.to_pandas()
        df = df.rename(columns=namemap)[namemap.values()]
        df['id_kic'] = df.id_kic.astype(int)
        df = add_prefix(df,'huber_')

    elif table=='huber14+cks':
        df = load_table('huber14',cache=1)
        cks = load_table('cks')
        df = pd.merge(df,cks)

    elif table=="huber14+cdpp":
        df = pd.read_hdf('data/kic.hdf')
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

        huber14 = load_table('huber14')
        stellar17 = cksphys.io.load_table('stellar17')
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

        return df

    elif table=='bruntt12':
        t = Table.read(
            'data/bruntt12/table3.dat',
            readme='data/bruntt12/ReadMe',format='ascii.cds'
        )
        namemap={'KIC':'id_kic','Teff':'steff','logg':'slogg','[Fe/H]':'smet'}
        df = t.to_pandas().rename(columns=namemap)[namemap.values()]

        df['id_kic'] = df['id_kic'].astype(int)
        df = add_prefix(df,'bruntt_')
    elif table=='bruntt12+cks':
        df = load_table('bruntt12')
        cks = load_table('cks')
        df = pd.merge(df,cks)
    elif table=='everett13':
        t = Table.read(
            'data/everett13/table2.dat',
            readme='data/everett13/ReadMe',
            format='ascii.cds'
            )
        namemap={'KIC':'id_kic','Teff':'steff','log(g)':'slogg','[Fe/H]':'smet'}
        df = t.to_pandas().rename(columns=namemap)[namemap.values()]
        df['id_kic'] = df['id_kic'].astype(int)
        df = add_prefix(df,'everett_')
    elif table=='everett13+cks':
        df = load_table('everett13')
        cks = load_table('cks')
        df = pd.merge(df,cks)

    elif table=='endl16':
        df = pd.read_csv('data/endl16_table2.txt',sep='\s') 
        namemap = {'KOI':'id_koi','Teff':'steff','logg':'slogg','FeH':'smet'}
        df = df.rename(columns=namemap)[namemap.values()]
        df = add_prefix(df,'endl_')
        df = df.convert_objects(convert_numeric=True)
    elif table=='endl16+cks':
        df = load_table('endl16')
        cks = load_table('cks')
        df = pd.merge(df,cks)

    elif table=='sm':
        df = pd.read_csv('data/specmatch_results.csv')
        df = df.groupby('name',as_index=False).last()
        df = df[df.name.str.contains('KIC|CK\d{5}|K\d{5}')]
        df['id_koi'] = None
        df['id_kic'] = None
        idx = df[df.name.str.contains('^K\d{5}$')].index

        df.ix[idx,'id_koi'] = df.ix[idx].name.str.slice(start=1).astype(int)

        idx = df[df.name.str.contains('^CK\d{5}$')].index
        df.ix[idx,'id_koi'] = df.ix[idx].name.str.slice(start=2).astype(int)

        idx = df[df.name.str.contains('^KIC')].index
        df.ix[idx,'id_kic'] = df.ix[idx].name.str.slice(start=3).astype(int)

        namemap = {'name':'id_name','id_kic':'id_kic','id_koi':'id_koi','teff':'steff','logg':'slogg','fe':'smet'}
        df = df.rename(columns=namemap)[namemap.values()]
        df = add_prefix(df,'sm_')
        df.index = df.id_name
        df = df.drop(['KIC3427720']) # Dropping because binary
        binaries = 'KIC9025370 KIC12317678'.split()
        df = df.drop(binaries)

    elif table=='sm+bruntt12':
        br = load_table('bruntt12')
        sm = load_table('sm')
        sm = sm.dropna(subset=['id_kic'])
        sm['id_kic'] = sm.id_kic.astype(int)
        br = br.dropna(subset=['id_kic'])
        br['id_kic'] = br.id_kic.astype(int)
        df = pd.merge(br,sm,on='id_kic')
    elif table=='sm+platinum':
        br = load_table('bruntt12')
        sm = load_table('sm')
        sm = sm.dropna(subset=['id_kic'])
        sm['id_kic'] = sm.id_kic.astype(int)
        br = br.dropna(subset=['id_kic'])
        br['id_kic'] = br.id_kic.astype(int)
        df = pd.merge(br,sm,on='id_kic')
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






