import pandas as pd
from astropy import constants as c 
from astropy import units as u
import numpy as np
from collections import Iterable
FLUX_EARTH = (c.L_sun / (4.0 * np.pi * c.au**2)).cgs.value
import os
import cksphys.io

from astropy.io import fits 
from cpsutils.pdplus import LittleEndian
import cksmet.cuts

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


    elif table=='nea':
        df = pd.read_csv('data/nea_2017-01-28.csv',comment='#')


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
        cuttypes = 'none faint diluted grazing longper badprad'.split()
        df = cksmet.cuts.add_cuts(df,cuttypes)
        
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


def radius_to_mass(radius):
    """
    Implement Weiss-Marcy Mass radius relationship
    """

    flux = 100 * FLUX_EARTH
    
    if isinstance(radius,Iterable):
        mass = map(radius_to_mass,radius)
        mass = np.array(mass)
        return mass
    
    if radius < 1.5: 
        rho = 2.43 + 3.39 * radius # g/cc WM14
        mass = (rho / 5.51) * radius**3
    elif 1.5 <= radius and radius <= 4.0: 
        mass = 2.69 * radius**0.93
    elif 4.0 < radius and radius < 9.0: 
        mass = 0.298 * radius**1.89 * flux**0.057
    elif 9.0 <= radius:
        mass = np.nan
    else:
        mass = np.nan
    return mass


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

def add_delta(pairs):
    """
    Compute the seperation between adjacent planet pairs interms of
    mutual hill-radii.

    """
    
    pairs['inner_mass'] = radius_to_mass(pairs.inner_iso_prad)
    pairs['outer_mass'] = radius_to_mass(pairs.outer_iso_prad)
    _delta = delta(
        pairs.inner_mass,
        pairs.outer_mass, 
        pairs.mass, 
        pairs.inner_smax, 
        pairs.outer_smax
    )
    pairs['delta'] = _delta
    return pairs

def hill_radius(mass1, mass2, stellar_mass, a1, a2):
    """
    mass1 (float): mass of inner planet (Earth Masses)
    mass2 (float): mass of outer planet (Earth Masses)
    stellar_mass (float): stellar mass (Solar Masses)
    a1 (float): semi-major axis (AU)


    """
    mass1 = np.array(mass1) * u.M_earth
    mass2 = np.array(mass2) * u.M_earth

    mass1 = mass1.to(u.M_sun).value
    mass2 = mass2.to(u.M_sun).value

    _hill_radius = ( 
        ((mass1 + mass2)/(3*stellar_mass))**(1.0/3.0) * 
        0.5*(a1 + a2)
    )
    return _hill_radius

def delta(mass1, mass2, stellar_mass, a1, a2):
    _hill_radius = hill_radius(mass1, mass2, stellar_mass, a1, a2)
    _delta = (a2 - a1)/_hill_radius
    return _delta




