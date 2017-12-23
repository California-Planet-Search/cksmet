import glob

import numpy as np
import pandas as pd

import cksmet.io
import cksmet.cuts

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
    lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1)
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
    lamo = cksmet.io.load_table('lamost-cal-cuts',cache=1)
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

def population(stub=False):
    lines = []
    df = cksmet.io.load_table('per-prad-population',cache=1)
    if stub:
        df = df.iloc[:10]
    df['i'] = range(len(df))
    for i, row in df.iterrows():
        line = r"{i:.0f} & {per:.3f} & {prad:.3f}  \\".format(**row)
        lines.append(line)

    return lines

