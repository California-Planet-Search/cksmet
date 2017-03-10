import cksmet.io
import cksmet.tables
import numpy as np
import cksphys
from cksmet.cuts import lamo_cuttypes, plnt_cuttypes, get_cut

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
    lamo = cksmet.io.load_table('lamost-dr2',cache=1)
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
    lamo = cksmet.io.load_table('lamost-dr2-cuts',cache=0)
    lamo = lamo[~lamo.isany]

    cks = cksmet.io.load_table('cks-cuts',cache=1)
    cks = cks.query('isany==False')
    
    quantiles = [0.01, 0.1,0.5,0.9,0.99]
    lamo = lamo.lamo_smet.quantile(quantiles)
    cks = cks.cks_smet.quantile(quantiles)

    for q, smet in lamo.iteritems():
        print r"{:.0f} & {:.2f} & {:.2f} \\".format(q*100, cks.ix[q], lamo.ix[q])


