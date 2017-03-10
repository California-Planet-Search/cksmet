"""
Module for dealing with cuts.
"""
import numpy as np
import sys
import inspect
import cksphys.io

cuttypes = 'none faint longper allfp badprad diluted grazing'.split()

class CutNotDwarf(object):
    cuttype = 'notdwarf'
    notstr = '$\log g = 3.9-5.0$ dex, $\mathregular{T}_\mathregular{eff} = 4700-6500$ (K)'
    def __call__(self,df):
        b1 = ~df.cks_slogg.between(3.9,5.0) | ~df.cks_steff.between(4700,6500)
        return b1

class CutAllFP(object):
    """Remove stars with a FP designation from a number of catalogs
    """
    cuttype = 'allfp'
    notstr = 'Not FP'
    texstr = 'Not a false positive'
    def cks(self,df):
        b1 = df['koi_disposition'].str.contains('FALSE POSITIVE')
        b2 = df['cks_fp']==1
        return b1 | b2


class CutNotCand(object):
    cuttype = 'notcand'
    notstr = 'Not FP (NEA)'
    texstr = 'Pl. Cand. (NEA)'
    def cks(self,df):
        return ~df['koi_disposition'].str.contains('CONFIRMED|CANDIDATE')

class CutFaint(object):
    cuttype = 'faint'
    notstr = '$Kp$ < 14.2'
    texstr = '$Kp$ < 14.2'
    def cks(self,df):
        return df['koi_kepmag'] > 14.2

class CutDiluted(object):
    cuttype= 'diluted'
    notstr = 'dilution < 5%'
    texstr = 'Radius correction factor < 5\%'
    def cks(self,df):
        return (df.furlan_rcorr_avg - 1) > 0.05

class CutGrazing(object):
    cuttype = 'grazing'
    notstr = '$b$ < 0.7' 
    texstr = '$b$ < 0.7'
    def cks(self,df):
        return df['koi_impact'] > 0.7 

class CutSubgiant(object):
    cuttype = 'subgiant'
    notstr = '$\log g < 3.9$'
    def cks(self,df):
        return df['cks_slogg'] < 3.9 

class CutLongPer(object):
    cuttype = 'longper'
    notstr = '$P$ < 350 d'
    texstr = '$P$ < 350 d'
    def cks(self,df):
        return df['koi_period'] > 350

class CutBadPrad(object):
    cuttype = 'badprad'
    notstr = '$\sigma(R_p) / R_p < 12\%$'
    texstr = 'Planet radius precision < 12\%'
    def cks(self,df):
        return df['iso_prad_err1']  / df['iso_prad'] > 0.12

class CutNone(object):
    cuttype = 'none'
    notstr = 'Full sample'
    texstr = 'Full sample'
    def cks(self,df):
        return np.zeros(len(df)).astype(bool)

## LAMOST Cuts
class CutLAMOSTNotDwarf(object):
    cuttype = 'lamonotdwarf'
    notstr = '$\log g = 3.9-5.0$ dex, $\mathregular{T}_\mathregular{eff} = 4700-6500$ (K)'
    def cks(self,df):
        b1 = ~df.lamo_slogg.between(3.9,5.0) | ~df.lamo_steff.between(4700,6500)
        return b1

clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)

def get_cut(cuttype):
    for name, obj in clsmembers:
        cut = obj()
        if getattr(cut,'cuttype')==cuttype:
            return cut
    assert False, "{} not defined".format(cuttype)

def add_cuts(df, cuttypes):
    isany = np.zeros(len(df))
    for cuttype in cuttypes:
        cut = get_cut(cuttype)
        key = 'is'+cuttype
        df[key]  = cut(df)
        isany += df[key].astype(int)

    df['isany'] = isany > 0 
    return df 

def count_cuts(df,cuttypes):
    """
    Return number of points that pass each cut
    """
    count = np.zeros(len(df))
    for cuttype in cuttypes:
        count += df['is'+cuttype].astype(int)
        print np.sum(count==0)

def table_of_cuts():
    """
    Apply cuts in sucession, count number of stars that pass
    """
    df =  cksphys.io.load_table(
        'cks+nea+iso-floor+huber-phot+furlan',cache=1,
        cachefn='../CKS-Physical/load_table_cache.hdf'
    )
    count = np.zeros(len(df))
    for cuttype in cuttypes:
        cut = get_cut(cuttype)
        npass = (cut.cks(df).astype(int)==False).sum()
        count += cut.cks(df).astype(int)
        npassall = np.sum(count==0)
        print r"{: <40} & {: <5} & {: <5} \\".format(cut.texstr,npass,npassall)

