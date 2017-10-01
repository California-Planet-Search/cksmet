import cksmet.io
import cksmet.tables
import numpy as np
import cksmet.cuts
from collections import OrderedDict
import pandas as pd

def cuts_planets():
    """
    Apply cuts in sucession, count number of stars that pass
    """
    cks =  cksmet.io.load_table('cks-cuts')
    cuttypes = cksmet.cuts.plnt_cuttypes
    lines = cut_statistics(cks, 'cks', cuttypes)
    return lines 

def cuts_field():
    """
    Apply cuts in sucession, count number of stars that pass
    """
    stars = cksmet.io.load_table('field-cuts',cache=1)
    cuttypes = cksmet.cuts.field_cuttypes
    lines = cut_statistics(stars, 'field', cuttypes)
    return lines

def cuts_lamost():
    """
    Apply cuts in sucession, count number of stars that pass
    """
    lamo = cksmet.io.load_table('lamost-dr2-cal-cuts',cache=1)
    cuttypes = cksmet.cuts.lamo_cuttypes
    lines = cut_statistics(lamo, 'lamo', cuttypes)
    return lines

def cut_statistics(df, sample, cuttypes):
    lines = []
    nall = len(df)
    count = np.zeros(nall)
    npassall = np.sum(count==0)
    npassallprev = np.sum(count==0)
    for cuttype in cuttypes:
        obj = cksmet.cuts.get_cut(cuttype)
        cut = obj(df,sample)
        b = cut.cut()
        npass = (b.astype(int)==False).sum()
        count += b.astype(int)
        npassall = np.sum(count==0)

        # pass all cuts
        f_pass_current = float(npassall) / float(npassallprev) 
        f_pass_all = float(npassall) / nall # pass all cuts
        line = r"{: <35} & {: <5} & {: <5} & {:.3f} \\".format(
            cut.texstr, npass, npassall, f_pass_current
        )
        lines.append(line)
        npassallprev = npassall

    return lines

def smet_stats():
    """
    Apply cuts in sucession, count number of stars that pass
    """
    lamo = cksmet.io.load_table('lamost-dr2-cal-cuts',cache=1)
    lamo = lamo[~lamo.isany]

    field = cksmet.io.load_table('field-cuts',cache=1)
    field = field.query('isany==False')

    cks = cksmet.io.load_table('cks-cuts',cache=1)
    cks = cks.query('isany==False')
    
    quantiles = [0.25, 0.5, 0.75]
    cks = cks.cks_smet
    lamo = lamo.lamo_smet
    field = field.m17_smet

    lamoq = lamo.quantile(quantiles)
    cksq = cks.quantile(quantiles)
    fieldq = field.quantile(quantiles)

    lines = []
    lines.append(r"Mean & {:.3f} & {:.3f} & {:.3f} \\".format(cks.mean(), field.mean(), lamo.mean()))
    lines.append(r"St. Dev. & {:.3f} & {:.3f}  & {:.3f} \\".format(cks.std(), field.std(), lamo.std()))
    func = lambda x : x.std()/np.sqrt(len(x))
    lines.append(r"St. Err. Mean & {:.3f} & {:.3f} & {:.3f} \\".format(
        func(cks), func(field), func(lamo)
    )
    )

    for q, smet in lamoq.iteritems():
        lines.append(r"{:.0f}\% & {:.3f} & {:.3f} & {:.3f} \\".format(q*100, cksq.ix[q], fieldq.ix[q], lamoq.ix[q]))
    
    return lines
def val_fit():
    fits = [
        'fit_per-sub-se',
        'fit_per-sup-se',
        'fit_per-sub-sn',
        'fit_per-sup-sn',

        'fit_smet-hot-se',
        'fit_smet-warm-se',

        'fit_smet-hot-sn',
        'fit_smet-warm-sn',

        'fit_smet-warm-ss',
        'fit_smet-hot-jup',
    ]
    lines = []
    for key in fits:
        _, prefix = key.split('_')
        fit = cksmet.io.load_object(key,cache=1)
        lines += fit.to_string(prefix=prefix+'-')
    return lines

def val_samp():
    d = OrderedDict()
    cache = 2
    cand = cksmet.io.load_table('cks-cuts',cache=cache)
    stars = cand.groupby('id_starname',as_index=False).first()

    d['n-cand-cks'] = "{}".format( len(cand) )
    d['n-stars-cks'] = "{}".format( len(stars) )

    cand = cand[~cand.isany]
    stars = cand.groupby('id_starname',as_index=False).first()

    d['n-cand-cks-pass'] = "{}".format( len(cand))
    d['n-stars-cks-pass'] = "{}".format( len(stars))
    d['n-stars-cks-pass-smet<-0.4'] = "{}".format( len(stars.query('cks_smet < -0.4')))

    field = cksmet.io.load_table('field-cuts',cache=1)
    d['n-stars-field'] = "{}".format( len(field))
    field = field[~field.isany]
    d['n-stars-field-pass'] = "{}".format( len(field))

    d['n-stars-field-pass-eff'] = cksmet.cuts.n_stars_field_pass_eff
    
    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    return lines


