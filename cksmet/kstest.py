import cksmet.io
from matplotlib.pylab import *
import pandas as pd
import seaborn as sns
import scipy.stats


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

    boxes = [
        dict(name="Hot Jupiters", per1=1.0, per2=10, prad1=8, prad2=24),
        dict(name="Warm Jupiters", per1=10, per2=100, prad1=8, prad2=24),
        dict(name="Cool Jupiters", per1=100, per2=350, prad1=8, prad2=24),

        dict(name="Hot Sub-Saturns", per1=1, per2=10, prad1=4, prad2=8),
        dict(name="Warm Sub-Saturns", per1=10, per2=100, prad1=4, prad2=8),
        dict(name="Cool Sub-Saturns", per1=100, per2=350, prad1=4, prad2=8),

        dict(name="Hot Sub-Neptunes", per1=1, per2=10, prad1=1.7, prad2=4),
        dict(name="Warm Sub-Neptunes", per1=10, per2=100, prad1=1.7, prad2=4),
        dict(name="Cool Sub-Neptunes", per1=100, per2=350, prad1=1.7, prad2=4),

        dict(name="Hot Super-Earths", per1=1, per2=10, prad1=1.0, prad2=1.7),
        dict(name="Warm Super-Earths", per1=10, per2=100, prad1=1.0, prad2=1.7),
        dict(name="Cool Super-Earths", per1=100, per2=350, prad1=1.0, prad2=1.7),

        dict(name="All Jupiters", per1=.1, per2=350, prad1=8, prad2=24),
        dict(name="All Sub-Saturns", per1=.1, per2=350, prad1=4, prad2=8),
        dict(name="All Sub-Neptunes", per1=.1, per2=350, prad1=1.7, prad2=4),
        dict(name="All Super-Earths", per1=.1, per2=350, prad1=1.0, prad2=1.7),
        dict(name="All Planets", per1=.1, per2=350, prad1=0.5, prad2=24),
    ]
    boxes = pd.DataFrame(boxes)
    cks = cksmet.io.load_table('cks-cuts',cache=1)
    cks = cks[~cks.isany]
    lamo = cksmet.io.load_table('lamost-dr2-cal-cuts',cache=2)
    lamo = lamo[~lamo.isany]
    #lamo = lamo.query('-1 < lamo_smet < 0.5')

    boxes = boxes.copy()
    #boxes = boxes.sort_values(by=['per1','prad1'])
    boxes2 = []
    

    lines = []
    for i,row in boxes.iterrows():
        cut = cks[
            cks.koi_period.between(row.per1,row.per2) &
            cks.iso_prad.between(row.prad1,row.prad2)
        ]

        d = calculate_statistics(cut.cks_smet,lamo.lamo_smet)
        d = dict(d, **row)
        lines.append( to_string(d))

        boxes2+=[d]

    boxes2 = pd.DataFrame(boxes2)
    #
    d = calculate_statistics(lamo.lamo_smet, lamo.lamo_smet)
    d['name'] = "LAMOST"
    s = "{name:s} & "
    s+="\\nodata & "
    s+="\\nodata & "
    s+="{n:.0f} & {mean:+.3f} & "
    s+="{sem:.3f} & "
    s+="\\nodata & \\nodata & \\nodata \\\\"
    s=s.format(**d)
    lines.append(s)
    return lines

def calculate_statistics(smet_plnt, smet_star):
    d = {}
    d['n'] = len(smet_plnt)
    d['mean'] = mean(smet_plnt)
    d['sem'] = std(smet_plnt) / np.sqrt(d['n'])
    stat, pval = scipy.stats.ks_2samp(smet_plnt, smet_star)
    d['kstest_pval'] = pval

    smet_star_mean = -0.005
    smet_sys = 0.02
    thresh = 0.01

    stat, pval = scipy.stats.ttest_1samp(smet_plnt, smet_star_mean)
    d['ttest_pval'] = pval
    d['ttest_stat'] = stat

    stat, pval = scipy.stats.ttest_1samp(smet_plnt, smet_star_mean-smet_sys)
    d['ttest_lo_pval'] = pval
    d['ttest_lo_stat'] = stat

    stat, pval = scipy.stats.ttest_1samp(smet_plnt, smet_star_mean+smet_sys)
    d['ttest_hi_pval'] = pval
    d['ttest_hi_stat'] = stat

    d = pd.Series(d)
    cols = ['ttest_pval','ttest_lo_pval','ttest_hi_pval']
    d['sig'] = np.sum(d[cols] < thresh)==len(cols)
    d['ttest_pval_max'] = d[cols].max()
    if d['sig']:
        d['sigstr'] = 'Y'
    else:
        d['sigstr'] = 'N'

    return d

def to_string(d):
    s = "{name:s} & "

#    if d['name'].count('All')==0:
    if True:
        s+="{per1:.0f}--{per2:.0f} & "
    else:
        s+="\\nodata & "
        
    s+="{prad1:.1f}--{prad2:.1f} & "
    s+="{n:.0f} & {mean:+.3f} & "
    s+="{sem:.3f} & "
    s+="\\num{{{ttest_stat:.2e}}} & < \\num{{{ttest_pval_max:.2e}}} & "
    s+="{sigstr:s} "
    s+="\\\\"
    s=s.format(**d)
    return s
