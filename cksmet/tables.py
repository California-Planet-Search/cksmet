import cksmet.io
import cksmet.tables
import numpy as np
import cksmet.cuts
from collections import OrderedDict
import pandas as pd
import numpy as np
import cksmet.io

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
        npassall_star = len(cut.df[count==0].id_kic.drop_duplicates())

        # pass all cuts
        f_pass_current = float(npassall) / float(npassallprev) 
        f_pass_all = float(npassall) / nall # pass all cuts
        line = r"{: <35} & {: <5} & {: <5} & {:.3f} & {: <5} \\".format(
            cut.texstr, npass, npassall, f_pass_current, npassall_star
        )
        lines.append(line)
        npassallprev = npassall

    return lines

def occurrence(stub=False):
    key = 'occur-per=0.25-prad=twoperoctave'
    occ = cksmet.io.load_object(key,cache=1)
    df = occ.df
    df = df.query('1 < perc < 300 and 0.5 < pradc < 32') 
    if stub:
        df = df.query('1 < pradc < 1.4') 

    df = df.sort_values(by=['pradc','perc'],ascending=[False,True])

    lines = []
    for i,row in df.iterrows():
        for k in 'rate_ul rate rate_err1 rate_err2'.split():
            row['p'+k] = row[k]*100

        if np.isnan(row.rate_err1):
            row['prate_tex'] = "< %(prate_ul).2f" % row
        else:
            row['prate_tex'] = "%(prate).2f^{+%(prate_err1).2f}_{%(prate_err2).2f}" % row
        if row.prob_det_mean < 0.25:
            row['prate_tex'] = r"\nodata"

        s=r""
        s+="{prad1:.2f} & {prad2:.2f} & "
        s+="{per1:.2f} & {per2:.2f} & "
        s+="{nplnt} & {prob_det_mean:.2f} & {ntrial:.1f} & "
        s+=r"${prate_tex}$ \\"
        s = s.format(**row)
        print s
        lines.append(s)

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
    lines.append(r"RMS & {:.3f} & {:.3f}  & {:.3f} \\".format(cks.std(), field.std(), lamo.std()))
    func = lambda x : x.std()/np.sqrt(len(x))
    lines.append(r"SEM & {:.3f} & {:.3f} & {:.3f} \\".format(
        func(cks), func(field), func(lamo)
    )
    )

    for q, smet in lamoq.iteritems():
        lines.append(r"{:.0f}\% & {:.3f} & {:.3f} & {:.3f} \\".format(q*100, cksq.ix[q], fieldq.ix[q], lamoq.ix[q]))
    
    return lines

def val_fit():
    fits = [
        'fit_per-all-se',
        'fit_per-all-sn',
        'fit_per-sub-se',
        'fit_per-sup-se',
        'fit_per-sub-sn',
        'fit_per-sup-sn',

        'fit_persmet-hot-se',
        'fit_persmet-hot-sn',
        'fit_persmet-hot-ss',
        'fit_persmet-hot-jup',

        'fit_persmet-warm-se',
        'fit_persmet-warm-sn',
        'fit_persmet-warm-ss',
        'fit_persmet-warm-jup',
    ]
    lines = []
    for key in fits:
        _, prefix = key.split('_')
        fit = cksmet.io.load_object(key,cache=1)
        lines += fit.to_string(prefix=prefix+'-')
    return lines

def val_samp(return_dict=False):
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

    # LAMOST Quantiles
    lamo = cksmet.io.load_table('lamost-dr2-cal-cuts',cache=1)
    lamo = lamo[~lamo.isany]
    quantiles = [0.25, 0.5, 0.75]
    lamoq = lamo.lamo_smet.quantile(quantiles)
    for q, smet in lamoq.iteritems():
        d['lamo-smet-{:.0f}'.format(q*100)] = "{:+.3f}".format(lamoq.ix[q])

    
    boxes = cksmet.kstest._boxes()
    boxes = pd.DataFrame(boxes)
    cks = cksmet.io.load_table('cks-cuts',cache=1)
    cks = cks[~cks.isany]

    boxes = boxes.copy()
    boxes2 = []
    lines = []
    for i,row in boxes.iterrows():
        cut = cks[
            cks.koi_period.between(row.per1,row.per2) &
            cks.iso_prad.between(row.prad1,row.prad2)
        ]
        _d = cksmet.kstest.calculate_statistics(cut.cks_smet,lamo.lamo_smet)
        _d = dict(_d, **row)
        d['{name} mean'.format(**_d)] = "{:+.2f}".format(_d['mean'])
        d['{name} pval'.format(**_d)] = "\\num{{{ttest_pval_max:.0e}}}".format(**_d)
        d['{name} sem'.format(**_d)] =  "{:.2f}".format(_d['sem'])
        d['{name} n'.format(**_d)] = "{n:.0f}".format(**_d)

    key = 'occur-per=0.25-prad=twoperoctave'
    occ = cksmet.io.load_object(key,cache=1)

    cut = occ.df.query('1 < perc < 10 and 8 < pradc < 24')
    print cut['rate'].sum()
    stats = cksmet.stats.sum_cells(cut.ntrial,cut.nplnt)
    d['rate-hot-jup'] = r"{%.2f}^{+%.2f}_{%.2f}" % (1e2*stats['rate'], 1e2*stats['rate_err1'], 1e2*stats['rate_err2'])
    d['rate-hot-jup-simple'] = r"{%.2f}" % (1e2*cut.rate.sum())

    if return_dict:
        return d

    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    return lines


import cksmet.kstest
