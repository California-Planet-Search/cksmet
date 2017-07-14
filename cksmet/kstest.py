import cksmet.io
from matplotlib.pylab import *
import pandas as pd
import seaborn as sns
import scipy.stats
boxes = [
    dict(name="Hot Jupiters", per1=1.0, per2=10, prad1=8, prad2=24),
    dict(name="Hot Small Planets", per1=0.3, per2=10, prad1=1, prad2=4),
    dict(name="Warm Small Planets", per1=10, per2=300, prad1=1, prad2=4),
    dict(name="Warm Sub-Saturns", per1=10, per2=300, prad1=4, prad2=8),
    dict(name="Warm Jupiters", per1=10, per2=300, prad1=8, prad2=16),
    dict(name="Warm Sub-Neptunes", per1=10, per2=300, prad1=1.7, prad2=4.0),
    dict(name="Warm Earths", per1=10, per2=300, prad1=1, prad2=1.7),
    dict(name="Hot Sub-Neptunes", per1=0.3, per2=10, prad1=1.7, prad2=4.0),
    dict(name="Hot Earths", per1=0.3, per2=10, prad1=1, prad2=1.7),
]
boxes = pd.DataFrame(boxes)

def plot_region():
    cks = cksmet.io.load_table('cks-cuts',cache=1)
    cks = cks[~cks.isany]
    loglog()
    plot(cks.koi_period, cks.iso_prad,'.')
    ax = gca()
    for i,row in boxes.iterrows():
        import matplotlib.patches as patches
        rect_kw = dict(fc='none',ec='k',lw=2) 
        
        corner = (row.per1,row.prad1)
        height = row.prad2 - row.prad1
        width = row.per2 - row.per1
        rect = patches.Rectangle(corner, width, height, **rect_kw)
        ax.add_patch(rect)
        text(corner[0]*1.05, corner[1]*1.05, row['name'], weight='bold')
        

    
def kstest_region():
    """
    Loop over the different regions and compute the KS test.
    """

    cks = cksmet.io.load_table('cks-cuts',cache=1)
    cks = cks[~cks.isany]

    lamo = cksmet.io.load_table('lamost-dr2-cal-cuts',cache=2)
    lamo = lamo[~lamo.isany]
    for i,row in boxes.iterrows():
        cut = cks[
            cks.koi_period.between(row.per1,row.per2) &
            cks.iso_prad.between(row.prad1,row.prad2)
        ]
 
        stat, pval = scipy.stats.ks_2samp(cut.cks_smet,lamo.lamo_smet)
        d = dict(row)
        d['pval'] = pval
        d['nplanets'] = len(cut)
        d['mean_smet'] = mean(cut.cks_smet)

        print "{name:s} & {per1:}--{per2:} & {prad1:}--{prad2:} & {nplanets:} & {mean_smet:+.2f} & \\num{{{pval:.2g}}} \\\\".format(**d)

    # All planets
    cut = cks
    stat, pval = scipy.stats.ks_2samp(cut.cks_smet,lamo.lamo_smet)
    d = dict(row)
    d['name'] = 'All Planets'
    d['pval'] = pval
    d['nplanets'] = len(cut)
    d['mean_smet'] = mean(cut.cks_smet)
    print "{name:s} & \\nodata & \\nodata & {nplanets:} & {mean_smet:+.2f} & \\num{{{pval:.2g}}} \\\\".format(**d)

    #
    cut = lamo
    d['name'] = "LAMOST"
    d['nplanets'] = len(lamo)
    d['mean_smet'] = mean(lamo.lamo_smet)
    print "{name:s} & \\nodata & \\nodata & {nplanets:} & {mean_smet:+.2f} & \\nodata \\\\".format(**d)

    # All planets
