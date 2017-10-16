"""
Module for dealing with cuts.
"""
import numpy as np
import sys
import inspect
import cksmet.io

# Effective number of field stars from which the planet sample was drawn.
n_stars_field_pass_eff = 33020

texdict ={
    'teff':'\\teff',
    'logg':'\\logg'
}

plotdict ={
    'teff':'\mathregular{T}_\mathregular{eff}',
    'logg':'\log g'
}

plnt_cuttypes = 'none faint badteffphot badloggphot longper allfp diluted grazing'.split()
lamo_cuttypes = 'none faint badteffphot badloggphot'.split()
field_cuttypes = 'none faint badteffphot badloggphot'.split()

samples = 'cks lamo field'.split()

class CutBase(object):
    cuttype = 'base'
    def __init__(self, df, sample):
        assert samples.count(sample)==1,'{} not supported sample'.format(sample)
        self.df = df
        self.sample = sample 
        self.nrows = len(df)

    def allpass(self):
        """A noop that lets all values pass
        """
        return np.zeros(self.nrows).astype(bool)

class CutNone(CutBase):
    cuttype = 'none'
    plotstr = 'Full sample'
    texstr = 'Full sample'
    def cut(self):
        return self.allpass()

class CutFaint(CutBase):
    cuttype = 'faint'
    plotstr = '$Kp$ < 14.2 mag'
    texstr = '$Kp$ < 14.2 mag'
    def cut(self):
        if self.sample=='cks':
            kepmag = self.df['kic_kepmag']
        elif self.sample=='lamo':
            kepmag = self.df['kic_kepmag']
        elif self.sample=='field':
            kepmag = self.df['kepmag']

        b = kepmag > 14.2
        return b

class CutTeffSpec(CutBase):
    """Remove stars where Teff falls outside allowed range
    """
    cuttype = 'badteffspec'
    texstr = r'${teff:}$ = $4700-6500$ K'.format(**texdict)
    plotstr = r'${teff:}$ = 4700$-$6500 K'.format(**plotdict)
    def cut(self):
        if self.sample=='cks':
            teff = self.df['cks_steff']
        elif self.sample=='lamo':
            teff = self.df['lamo_steff']

        b = ~teff.between(4700,6500)
        return b

class CutTeffPhot(CutBase):
    """Remove stars where Teff falls outside allowed range
    """
    cuttype = 'badteffphot'
    texstr = r'${teff:}$ = $4700-6500$ K'.format(**texdict)
    plotstr = r'${teff:}$ = 4700$-$6500 K'.format(**plotdict)
    def cut(self):
        if self.sample=='cks':
            teff = self.df['m17_steff']
        elif self.sample=='lamo':
            teff = self.df['m17_steff']
        elif self.sample=='field':
            teff = self.df['m17_steff']

        b = ~teff.between(4500,6500)
        return b

class CutLoggSpec(CutBase):
    """Remove stars where logg falls outside allowed range
    """
    cuttype = 'badloggspec'
    texstr = '${logg:}$ = $3.9-5.0$ dex'.format(**texdict)
    plotstr = 'log g = 3.9$-$5.0 dex'.format(**plotdict)
    def cut(self):
        if self.sample=='cks':
            logg = self.df['cks_slogg']
        elif self.sample=='lamo':
            logg = self.df['lamo_slogg']
        elif self.sample=='field':
            logg = self.df['huber_slogg']
        b = ~logg.between(3.9,5.0)
        return b

class CutLoggPhot(CutBase):
    """Remove stars where logg falls outside allowed range
    """
    cuttype = 'badloggphot'
    texstr = '${logg}$ = $3.9-5.0$ dex'.format(**texdict)
    plotstr = 'log g = 3.9$-$5.0 dex'.format(**plotdict)
    def cut(self):
        if self.sample=='cks':
            logg = self.df['m17_slogg']
        elif self.sample=='lamo':
            logg = self.df['m17_slogg']
        elif self.sample=='field':
            logg = self.df['m17_slogg']
        b = ~logg.between(3.9,5.0)
        return b

class CutAllFP(CutBase):
    """Remove stars with a FP designation from a number of catalogs
    """
    cuttype = 'allfp'
    plotstr = 'Not FP'
    texstr = 'Not a false positive'
    def cut(self):
        if self.sample=='cks':
            b1 = self.df['koi_disposition'].str.contains('FALSE POSITIVE')
            b2 = self.df['cks_fp']==1
            return b1 | b2
        elif self.sample=='field':
            return self.allpass()

class CutLowPPrad(CutBase):
    cuttype = 'lowpprad'
    texstr = '$\sigma(R_p) / R_p < 12\%$'
    plotstr = texstr
    def cut(self):
        if self.sample=='cks':
            return self.df['iso_prad_err1']  / self.df['iso_prad'] > 0.12
        elif self.sample=='field':
            return self.allpass()

class CutLowPSrad(CutBase):
    cuttype = 'lowpsrad'
    texstr = '$\sigma(R_{\star}) / R_{\star} < 20\%$'
    plotstr = texstr
    def cut(self):
        if self.sample=='cks':
            return self.df['iso_srad_err1']  / self.df['iso_srad'] > 0.20
        elif self.sample=='field':
            return self.allpass()

class CutDiluted(CutBase):
    cuttype= 'diluted'
    plotstr = 'dilution < 5%'
    texstr = 'Radius correction factor < 5\%'
    def cut(self):
        if self.sample=='cks':
            return (self.df.furlan_rcorr_avg - 1) > 0.05
        elif self.sample=='field':
            return self.allpass()

class CutGrazing(CutBase):
    cuttype = 'grazing'
    plotstr = 'Not grazing' 
    texstr = 'Not a grazing transit ($b$ < 0.9)'
    def cut(self):
        if self.sample=='cks':
            return self.df['koi_impact'] > 0.9 
        elif self.sample=='field':
            return self.allpass()

class CutLongPer(CutBase):
    cuttype = 'longper'
    plotstr = '$P$ < 350 days'
    texstr = '$P$ < 350 days'
    def cut(self):
        if self.sample=='cks':
            return self.df['koi_period'] > 350
        elif self.sample=='lamo':
            return self.allpass()
        elif self.sample=='field':
            return self.allpass()

class CutSmallPrad(CutBase):
    cuttype = 'smallprad'
    plotstr = '$R_P$ > 0.5 $R_{\oplus}$'
    texstr = plotstr
    def cut(self):
        if self.sample=='cks':
            b1 = self.df['iso_prad'] < 0.5
            b2 = self.df['iso_prad'].isnull()
            return b1 | b2
        elif self.sample=='lamo':
            return self.allpass()
        elif self.sample=='field':
            return self.allpass()

class CutNotDwarf(CutBase):
    """
    """
    cuttype = 'notdwarf'
    plotstr = 'Dwarf star'
    texstr = 'Dwarf star'
    def cut(self):
        if self.sample=='cks':
            srad = self.df['iso_srad']
            logsrad = np.log10(srad)
            steff = self.df['iso_steff']
            b1 = logsrad > 0.00025 * (steff - 5500) + 0.2
            #b2 = (srad < 0.8) | (srad > 1.2)
            b2 = (srad < 0.8) 
            return b1 | b2
        elif self.sample=='field':
            srad = self.df['huber_srad']
            logsrad = np.log10(srad)
            steff = self.df['huber_steff']
            b1 = logsrad > 0.00025 * (steff - 5500) + 0.2
            #b2 = (srad < 0.8) | (srad > 1.2)
            b2 = (srad < 0.8) 
            return b1 | b2

clsmembers = inspect.getmembers(sys.modules[__name__], inspect.isclass)

def get_cut(cuttype):
    for name, obj in clsmembers:
        if getattr(obj,'cuttype')==cuttype:
            return obj

    assert False, "{} not defined".format(cuttype)

def add_cuts(df, cuttypes, sample):
    isany = np.zeros(len(df))
    for cuttype in cuttypes:
        obj = get_cut(cuttype)
        cut = obj(df,sample)
        key = 'is'+cuttype
        df[key] = cut.cut()
        isany += df[key].astype(int)

    df['isany'] = isany > 0 
    return df 
