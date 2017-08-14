import cksmet.io
import cksmet.tables
import numpy as np
import cksphys
from cksmet.cuts import lamo_cuttypes, plnt_cuttypes, get_cut
from collections import OrderedDict
import cksphys.io
import pandas as pd
from collections import OrderedDict



def cuts_planets():
    """
    Apply cuts in sucession, count number of stars that pass
    """
    cks =  cksphys.io.load_table(
        'cks+nea+iso-floor+huber-phot+furlan',cache=1,
        cachefn='../CKS-Physical/load_table_cache.hdf'
    )
    print_cut_statistics(cks, 'cks', plnt_cuttypes)

def cuts_field():
    """
    Apply cuts in sucession, count number of stars that pass
    """
    stars = cksmet.io.load_table('huber14+cdpp',cache=1)
    print_cut_statistics(stars, 'field', plnt_cuttypes)

def cuts_lamost():
    """
    Apply cuts in sucession, count number of stars that pass
    """
    lamo = cksmet.io.load_table('lamost-dr2-cal',cache=1)
    print_cut_statistics(lamo, 'lamo', lamo_cuttypes)

def print_cut_statistics(df, sample, cuttypes):
    nall = len(df)
    count = np.zeros(nall)
    npassall = np.sum(count==0)
    npassallprev = np.sum(count==0)
    for cuttype in cuttypes:
        obj = get_cut(cuttype)
        cut = obj(df,sample)
        b = cut.cut()
        npass = (b.astype(int)==False).sum()
        count += b.astype(int)
        npassall = np.sum(count==0)

        # pass all cuts
        f_pass_current = float(npassall) / float(npassallprev) 
        f_pass_all = float(npassall) / nall # pass all cuts
        print r"{: <35} & {: <5} & {: <5} & {:.2f} \\".format(
            cut.texstr, npass, npassall, f_pass_current
        )
        npassallprev = npassall

def smet_dist_lamost():
    """
    Apply cuts in sucession, count number of stars that pass
    """
    lamo = cksmet.io.load_table('lamost-dr2-cal-cuts',cache=0)
    lamo = lamo[~lamo.isany]

    cks = cksmet.io.load_table('cks-cuts',cache=1)
    cks = cks.query('isany==False')
    
    quantiles = [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]
    cks = cks.cks_smet
    lamo = lamo.lamo_smet
    lamoq = lamo.quantile(quantiles)
    cksq = cks.quantile(quantiles)

    print r"% statistic, CKS, LAMOST"
    print r"Mean & {:.3f} & {:.3f} \\".format(cks.mean(), lamo.mean())
    print r"St. Dev. & {:.3f} & {:.3f} \\".format(cks.std(), lamo.std())
    print r"St. Err. Mean & {:.3f} & {:.3f} \\".format(
        cks.std()/np.sqrt(len(cks)), lamo.std()/np.sqrt(len(lamo))
    )


    for q, smet in lamoq.iteritems():
        print r"{:.0f}\% & {:.3f} & {:.3f} \\".format(q*100, cksq.ix[q], lamoq.ix[q])


def print_fit_stats():
    fits = [
        'fit_per-sub-se',
        'fit_per-sup-se',
        'fit_per-sub-sn',
        'fit_per-sup-sn',
        'fit_per-sup-ss',

        'fit_smet-hot-se',
        'fit_smet-warm-se',

        'fit_smet-hot-sn',
        'fit_smet-warm-sn',

        'fit_smet-warm-ss',
        'fit_smet-hot-jup',
    ]

        
    for key in fits:
        _, prefix = key.split('_')
        fit = cksmet.io.load_fit(key,cache=1)
        fit.print_parameters(prefix+'-') 


def stats():
    d = OrderedDict()
    cache = 1 
    cand = cksphys.io.load_table('cks-cuts',cache=cache)
    stars = cand.groupby('id_starname',as_index=False).first()

    d['n-cand-cks'] = "{}".format( len(cand) )
    d['n-stars-cks'] = "{}".format( len(stars) )

    cand = cand[~cand.isany]
    stars = cand.groupby('id_starname',as_index=False).first()

    d['n-cand-cks-pass'] = "{}".format( len(cand))
    d['n-stars-cks-pass'] = "{}".format( len(stars))

    for k, v in d.iteritems():
        print r"{{{}}}{{{}}}".format(k,v)



