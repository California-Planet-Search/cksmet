from collections import OrderedDict
import glob

import pandas as pd
import numpy as np

import cksmet.io
import cksmet.tables
import cksmet.ttest

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

    df = cksmet.io.load_table('lamost-cks-calibration-sample-noclip',cache=1)
    d['n-cks-lamo-tot'] = len(df)

    # LAMOST calibration
    df = cksmet.io.load_table('lamost-cks-calibration-sample',cache=1)
    d['n-cks-lamo-cal'] = len(df)
    d['n-cks-lamo-out'] = d['n-cks-lamo-tot'] - d['n-cks-lamo-cal']


    # LAMOST 
    lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1)
    d['n-stars-lamo'] = "{}".format( len(lamo) )
    lamo = lamo[~lamo.isany]
    d['n-stars-lamo-pass'] = "{}".format( len(lamo))

    quantiles = [0.25, 0.5, 0.75]
    lamoq = lamo.lamo_smet.quantile(quantiles)
    for q, smet in lamoq.iteritems():
        d['lamo-smet-{:.0f}'.format(q*100)] = "{:+.3f}".format(lamoq.ix[q])

    boxes = cksmet.ttest._boxes()
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
        _d = cksmet.ttest.calculate_statistics(cut.cks_smet,lamo.lamo_smet)
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

    x = dict(fe=0.0)
    dx = dict(fe=0.1)
    def hyperplane_coeff_samples(prefix, param, x,dx,):
        fL = glob.glob(prefix)
        samples = []
        for f in fL[:100]:
            cal = cksmet._calibrate.read_fits(f,param)
            out = cal.hyperplane_coeff(x,dx)
            out['param'] = param
            samples.append(out)
        return samples

    L1 = hyperplane_coeff_samples('bootstrap/L1*','fe',x,dx)
    L2 = hyperplane_coeff_samples('bootstrap/L2*','fe',x,dx)

    L1 = pd.DataFrame(L1)
    L2 = pd.DataFrame(L2)

    d['L2-c0'] = "{:.3f}".format(L2['c0'].mean())
    d['L2-c1'] = "{:.3f}".format(L2['cfe'].mean())
    d['L2-c0-err'] = "{:.3f}".format(L2['c0'].std())
    d['L2-c1-err'] = "{:.3f}".format(L2['cfe'].std())

    xeval = -0.5
    samp = L2.c0 + L2.cfe * (xeval - x['fe'])/dx['fe']
    d['L2-val-m0.5'] = "{:+.3f} \pm {:.3f}".format(samp.mean(),samp.std())

    xeval = +0.5
    samp = L2.c0 + L2.cfe * (xeval - x['fe'])/dx['fe']
    d['L2-val-p0.5'] = "{:+.3f} \pm {:.3f}".format(samp.mean(),samp.std())

    d['L1-c0'] = "{:.3f}".format(L1['c0'].mean())
    d['L1-c1'] = "{:.3f}".format(L1['cfe'].mean())
    d['L1-c0-err'] = "{:.3f}".format(L1['c0'].std())
    d['L1-c1-err'] = "{:.3f}".format(L1['cfe'].std())

    p, nplanets = cksmet.population.load_pinky()
    nstars = cksmet.population.nstars
    d['pop-nstars'] = nstars
    d['pop-nplanets'] = nplanets
    d['pop-intrate'] = "{:.1f}".format(100.0 * nplanets / nstars)

    lines = []
    for k, v in d.iteritems():
        line = r"{{{}}}{{{}}}".format(k,v)
        lines.append(line)

    if return_dict:
        return d
    return lines

