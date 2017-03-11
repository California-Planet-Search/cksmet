from matplotlib.pylab import *
import seaborn as sns
import pandas as pd

import cksphys.io
import cksmet.io
import cksmet.cuts
from collections import Iterable
from matplotlib.ticker import FuncFormatter, MaxNLocator

sns.set_style('whitegrid')
sns.set_color_codes()

texrp = '\mathregular{R}_\mathregular{P}' 
texre = '\mathregular{R}_\mathregular{E}' 


rpticks = [0.2, 0.3, 0.4, 0.5, 0.7, 1, 2, 3, 4, 5, 7, 10, 20]

from cksmet.plotting.samples import add_anchored

def prad_fe_label():
    # median planet radius error
    xlabel('Planet size (Earth-radii)')
    ylabel('[Fe/H]')
    xticks(rpticks,rpticks)

def prad_fe_errbar():
    cks = cksmet.io.load_table('cks-cuts')
    ebarx = 14
    ebary = -0.7
    xferr = (cks.iso_prad_err1 / cks.iso_prad).median()    
    xerr = [[ ebarx - ebarx/(1+xferr) ], [ebarx*(1+xferr) - ebarx]]
    yerr = [[0.04],[0.04]]
    errorbar([ebarx], [ebary], yerr=yerr, xerr=xerr,fmt='o')
    annotate( 
        xy=(ebarx,ebary), s='Median   \nUncert.   ',ha='right', 
        va='baseline'
    )

def bins_to_xerr(bin0,binc,bin1):
    xerr = [[np.array(binc - bin0)], [np.array(bin1 - binc)] ]
    xerr = np.vstack(xerr)
    return xerr

def prad_fe():
    """Plot of planet radius vs stellar metallicty"""
    cks = cksmet.io.load_table('cks-cuts')
    cks = cks[~cks.isany]
#    bins = [0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 5.7, 8.0, 11.3, 16]
    bins = [0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 8.0, 16]
#    bins = [0.5, 1.7, 4.0, 8.0, 16]
    cksbin = cksmet.io.table_bin(cks, bins)

    figure(figsize=(6,4))

    semilogx()
    plot(cks.iso_prad,cks.cks_smet,'.')
    prad_fe_label()
    prad_fe_errbar()

    i = 0 
    for _i, row in cksbin.iterrows():
        x = np.sqrt(row.bin0*row.bin1)
        xerr = [[x - row.bin0], [row.bin1 - x] ]
        yerr = row.fe_mean_err
        y = row.fe_mean
        errorbar(x,y,xerr=xerr,yerr=yerr,color='r',zorder=10)
        i+=1

    xlim(0.25,25)
    ylim(-0.9,0.7)


def prad_fe_percentiles():
    """Plot of planet radius vs stellar metallicty"""
    cks = cksmet.io.load_table('cks-cuts')
    cks = cks[~cks.isany]

    bins = [0.7, 1.0, 1.4, 2.0, 2.8, 4.0, 8.0, 16]
    #cks = cks.query(' -0.3 < cks_smet < 0.3')
    cksbin = cksmet.io.table_bin(cks, bins)

    figure(figsize=(6,4))
    semilogx()
    plot(cks.iso_prad,cks.cks_smet,'.',color='LightGray')
    prad_fe_label()

    i = 0 
    for _i, row in cksbin.iterrows():
        x = np.sqrt(row.bin0*row.bin1)
        xerr = [[x - row.bin0], [row.bin1 - x] ]
        yerr = row.fe_mean_err
        y = row.fe_mean
        errorbar(x, row.fe_mean, xerr=xerr,color='r',zorder=10)
        errorbar(x,row.fe_p01,xerr=xerr,color='b',zorder=10)
        errorbar(x,row.fe_p50,xerr=xerr,color='b',zorder=10)
        errorbar(x,row.fe_p99,xerr=xerr,color='b',zorder=10)
        i+=1

    xlim(0.25,25)
    ylim(-0.9,0.7)

def prad_fe_mean():
    """Plot of planet radius vs stellar metallicty"""
    semilogx()
    cksbin = cksmet.io.load_table('cksbin-fe')
    xerr =  np.vstack([[cksbin.binc - cksbin.bin0], [cksbin.bin1 - cksbin.binc] ])

    errorbar(cksbin.binc,cksbin.fe_mean,yerr=cksbin.fe_mean_err,xerr=xerr,fmt='o')
    errorbar(cksbin.binc,cksbin.fe_p50,yerr=cksbin.fe_mean_err,xerr=xerr,fmt='o')
    xlim()
    xticks(rpticks,rpticks)

def prad_fe_multiplicity():
    """Plot of planet radius vs stellar metallicty multiplicity labeled"""

    cks = cksmet.io.load_table('cks')
    figure(figsize=(8,6))
    semilogx()
    plot(cks.iso_prad,cks.cks_smet,'.')


def period_fe():
    """Plot of planet orbital period vs stellar metallicty 
    """
    figure(figsize=(8,6))
    semilogx()
    plot(cks.koi_period, cks.cks_smet,'.')
    cut = cks.query('nplanets > 5')
    plot(cut.koi_period, cut.cks_smet,'.r')
    #prad_fe_label()


def prad_fe_fit():
    cksbin = cksmet.io.load_table('cksbin-fe')
    x = log10(cksbin.binc)
    y = cksbin.fe_mean
    yerr = cksbin.fe_mean_err
    import lmfit

    p0 = lmfit.Parameters()
    p0.add('p1',0.1)
    p0.add('p0',0.1)

    def model_linear(params,x):
        return params['p1'] * x + params['p0']

    def model_step(params, x):
        if isinstance(x,Iterable):
            out = map(lambda x : model_step(params,x), x)
            out = np.array(out)
            return out

        if x < log10(4):
            return params['p0']
        elif x >=  log10(4):
            return params['p1']

    def resid(params, model, x, y, yerr):
        return (y - model(params,x))/yerr

    res_linear = lmfit.minimize(resid,p0,args=(model_linear, x,y,yerr,))
    p1_linear = res_linear.params
    print "linear model"
    lmfit.report_fit(res_linear)

    res_step = lmfit.minimize(resid,p0,args=(model_step, x,y,yerr,))
    p1_step = res_linear.params
    print "step model"
    lmfit.report_fit(res_step)


    print "BIC(linear) - BIC(step) = {}".format(res_linear.bic - res_step.bic)

    semilogx()
    xerr = bins_to_xerr(cksbin.bin0,cksbin.binc,cksbin.bin1)
    print xerr
    errorbar(10**x,y,xerr=xerr, yerr=yerr,fmt='o')
    plot(10**x,model_linear(p1_linear,x))

    prad_fe_label()
    xlim(0.5,16)


def smax_fe4():
    cks = cksmet.io.load_table('cks-cuts')
    cks = cks[~cks.isany]
    fig = figure(figsize=(8.5,11))
    fig, axL = subplots_4()
    sns.husl_palette(5,h=.7)
    pal = sns.color_palette("muted", 4)
    sns.set_palette(pal)

    radius_cuts = [
        (0.5, 1.8),
        (1.8, 4),
        (4.0, 8.0),
        (8.0, 16.0),
    ]

    sca(axL[0])
    semilogx()
    plot(cks.iso_sma,cks.cks_smet,'.')
    
    for i in range(4):
        sca(axL[i+1])
        radius_cut = radius_cuts[i]
        print radius_cut
        label = r"{0:.1f}-{1:.1f} $\mathregular{{R}}_\mathregular{{E}}$".format(
            *radius_cut
        )
        cut = cks[cks.iso_prad.between(*radius_cut)]
        semilogx()
        plot(cks.iso_sma,cks.cks_smet,'.',color='LightGray',label='')
        plot(cut.iso_sma,cut.cks_smet,'.',label=label)
        xt = [0.003,0.01,0.03,0.1,0.3,1,3]
        xticks(xt,xt)
        legend(loc='lower right',title='',frameon=True)

    setp(axL,xlim=(0.004,3),ylim=(-1.0,0.6))
    setp(axL[0],xlabel='a [AU]',ylabel='[Fe/H]')
    fig.set_tight_layout(True)


def subplots_4(*args, **kwargs):
    fig = gcf()
    shape = (4,2)
    ax1 = plt.subplot2grid(shape, (0, 0) , colspan=2, rowspan=2)
    ax2 = plt.subplot2grid(shape, (2, 0))
    ax3 = plt.subplot2grid(shape, (2, 1))
    ax4 = plt.subplot2grid(shape, (3, 0))
    ax5 = plt.subplot2grid(shape, (3, 1))
    axL = [ax1,ax2,ax3,ax4,ax5]
    return fig, axL

def labelfig(text):
    add_anchored(
        text, loc=2, prop=dict(fontweight='bold',size='medium'), 
        frameon=False
    )


def nplanets_fe():
    """Plot of planet radius vs stellar metallicty multiplicity labeled"""

    fig = figure(figsize=(4,4)) 
    shape = (3,1)
    ax1 = plt.subplot2grid(shape, (0, 0) , colspan=1, rowspan=2)
    ax2 = plt.subplot2grid(shape, (2, 0))
    axL = [ax1,ax2]
    sca(ax1)
    cks = cksmet.io.load_table('cks')
    plot(cks.nplanets, cks.cks_smet, '_', lw=1, mew=1, ms=2, alpha=0.5)
    xlim(0,8)
    cksbin = cksmet.io.load_table('cksbin-nplanets')
    print cksbin
    
    for ax in axL:
        sca(ax)
        errorbar(
            cksbin.multiplicity, cksbin.fe_mean, yerr=cksbin.fe_mean_err,
            fmt='_',mew=2,color='g'
        )
        
    fig.set_tight_layout(True)
    setp(axL,xlim=(0,8),ylabel='[Fe/H]')
    xlabel('Multiplicity')

def radius_enhancement(binkey='pl_rade', bins=[1,2,3,4,8]):
    df = cksmet.io.load_table('hadden+cks+nea',cache=1)
    query = (
    "pl_radeerr1 / pl_rade < 0.5 "
    + "and -pl_radeerr2 / pl_rade < 0.5 "
    + "and pl_masseerr1 / pl_masse < 0.5 "
    + "and -pl_masseerr2 / pl_masse < 0.5"
    + "and pl_masse < 100"
    )
    df = df.query(query)

    df['pl_logmasse'] = np.log10(df.pl_masse)
    df['pl_lograde'] = np.log10(df.pl_rade)

    febins = [-0.5,0.05,0.5]
    loglog()
    for i in range(len(febins)-1):
        fe1 = febins[i]
        fe2 = febins[i+1]

        cut = df[df.cks_smet.between(fe1,fe2)]
        if len(cut)==0:
            continue

        x = cut.pl_masse
        y = cut.pl_rade
        xerr = cut.pl_masseerr1
        yerr = cut.pl_radeerr1

        errorbar(
            x,y,xerr=xerr,yerr=yerr ,fmt='.',label='[{} -- {}]'.format(fe1,fe2)
        )
    xlabel('Planet Mass (Earth-masses)')
    ylabel('Planet Size (Earth-radii)')

    legend()
    cut = df.query('pl_logmasse < 2')
    p1 = np.polyfit(cut.pl_logmasse,cut.pl_lograde,1)
    p1 = [0.7,-0.2]

    pl_logmassei = linspace(0,3,100)
    pl_logradefit = np.polyval(p1,pl_logmassei)
    plot(10**pl_logmassei,10**pl_logradefit)

    fig, axL = subplots(nrows=2,ncols=3,sharex=True,sharey=True)
    axL = axL.flatten()

    df['pl_lograde_diff'] = df.pl_lograde - np.polyval(p1,df.pl_logmasse)

    semilogy()

    for i in range(len(bins)-1):
        sca(axL[i])
        

        bin1 = bins[i]
        bin2 = bins[i+1]
        cut = df[df[binkey].between(bin1,bin2)]

        label = binkey + " = {} -- {}".format(bin1,bin2)
        plot(cut.cks_smet,10**cut.pl_lograde_diff,'.',)
        title(label)
        yt = [0.3,0.5,1,2,3]
        yticks(yt,yt)
        xt = [-0.4,-0.2,0.0,0.2,0.4]
        xticks(xt,xt)

    sca(axL[3])
    xlabel('metalicity')
    ylabel('radius enhancement')
    xlim(-0.5,0.5,)
    ylim(0.2,5)


def mass_enhancement(binkey='pl_rade', bins=[1,2,3,4,8]):
    df = cksmet.io.load_table('hadden+cks+nea',cache=1)
    query = (
    "pl_radeerr1 / pl_rade < 0.5 "
    + "and -pl_radeerr2 / pl_rade < 0.5 "
    + "and pl_masseerr1 / pl_masse < 0.5 "
    + "and -pl_masseerr2 / pl_masse < 0.5"
    + "and pl_masse < 100"
    )
    df = df.query(query)

    df['pl_logmasse'] = np.log10(df.pl_masse)
    df['pl_lograde'] = np.log10(df.pl_rade)

    febins = [-0.5,0.05,0.5]
    loglog()
    for i in range(len(febins)-1):
        fe1 = febins[i]
        fe2 = febins[i+1]

        cut = df[df.cks_smet.between(fe1,fe2)]
        if len(cut)==0:
            continue

        yerr = cut.pl_masseerr1
        xerr = cut.pl_radeerr1
        y = cut.pl_masse
        x = cut.pl_rade

        errorbar(
            x,y,xerr=xerr,yerr=yerr ,fmt='.',label='[{} -- {}]'.format(fe1,fe2)
        )
    xlabel('Planet Radius (Earth-radii)')
    ylabel('Planet Mass (Earth-masses)')

    legend(loc='best')
    cut = df.query('pl_logmasse < 2')
    p1 = np.polyfit(cut.pl_lograde,cut.pl_logmasse,1)
    p1 = [1.3,+0.2]

    pl_logradei = linspace(0,3,100)
    pl_logmassefit = np.polyval(p1,pl_logradei)
    plot(10**pl_logradei,10**pl_logmassefit)
    xlim(0,20)
    ylim(0,100)
    xt = [1,2,3,4,5,7,10,20,]
    xticks(xt,xt)
    yt = [1,3,10,30,100]
    yticks(yt,yt)

    fig, axL = subplots(nrows=2,ncols=3,sharex=True,sharey=True)
    axL = axL.flatten()

    df['pl_logmasse_diff'] = df.pl_logmasse - np.polyval(p1,df.pl_lograde)

    semilogy()

    for i in range(len(bins)-1):
        sca(axL[i])
        bin1 = bins[i]
        bin2 = bins[i+1]
        cut = df[df[binkey].between(bin1,bin2)]

        label = binkey + " = {} -- {}".format(bin1,bin2)
        plot(cut.cks_smet,10**cut.pl_logmasse_diff,'.',)
        title(label)
        yt = [0.3,0.5,1,2,3]
        yticks(yt,yt)
        xt = [-0.4,-0.2,0.0,0.2,0.4]
        xticks(xt,xt)

    sca(axL[3])
    xlabel('metalicity')
    ylabel('mass enhancement')
    ylim(0.1,10)
    xlim(-0.5,0.5)

def prad_hist_stacked(prad_bins=[0.5,1.0,2.0,4.0,8.0,16]):
    lamo = cksmet.io.load_table('lamost-dr2-cuts',cache=1)
    cks = cksmet.io.load_table('cks-cuts',cache=1)

    lamo = lamo.query('isany==False')
    cks = cks.query('isany==False')

    smet_bins = arange(-0.8,0.81,0.1)
    nrows = len(prad_bins) - 1 
    ncols = 2
    width = 7.5
    height = 9
    fig, axL = subplots(
        nrows=nrows,ncols=ncols,figsize=(width,height),sharex=True
    )

    # Lamost taps out at kepmag = 14, so let's restrict
    i = 0 
    cut = cks
    for i in range(nrows):
        histogramkw = dict(bins=smet_bins,)
        histkw = dict(histtype='step',lw=2, bins=smet_bins)
        prad1 = prad_bins[i]
        prad2 = prad_bins[i+1]
        irow = nrows-1-i
        
        axhist = axL[irow,0]
        axcum = axL[irow,-1]

        sca(axhist)

        cut = cks[cks.iso_prad.between(prad1,prad2)]

        # LAMOST HIST
        counts,_ = histogram(lamo.lamo_smet,**histogramkw)
        total = len(lamo)
        weights = 1.0*counts / total
        hist(smet_bins[:-1], weights=weights, color='Gray',**histkw)

        # CKS HIST
        counts,_ = histogram(cut.cks_smet,**histogramkw)
        total = len(cut)
        weights = 1.0*counts / total
        hist(smet_bins[:-1], weights=weights, color='b',**histkw)

        sca(axcum)
        histkw = dict(
            bins=smet_bins, normed=True,cumulative=True,histtype='step',lw=2,zorder=10
        )
        hist(lamo.lamo_smet, color='Gray', **histkw)
        hist(cut.cks_smet, color='b', **histkw)

        sca(axhist)
        x1 = 0.02
        y1 = 0.98
        y2 = 0.72
        textkw = dict(transform=axhist.transAxes, fontsize='medium',va='top')
        text(x1, y1,'CKS\n${}$ = {}$-${} ${}$'.format(texrp,prad1,prad2,texre),color='b',**textkw)
        text(x1, y2,'LAMOST',color='Gray',**textkw)

        sca(axcum)
        textkw['transform'] = axcum.transAxes
        text(x1, y1,'CKS\n${}$ = {}$-${} ${}$'.format(texrp,prad1,prad2,texre),color='b',**textkw)
        text(x1, y2,'LAMOST',color='Gray',**textkw)


    setp(axL[:,0],xlim=(-1.0,0.6),ylim=(0,0.3))
    setp(axL[:,1],xlim=(-1.0,0.6),ylim=(0,1.0))
    setp(axL[-1,:],xlabel='[Fe/H]')

    sca(axL[2,0])
    ylabel('Fraction of Candidates per Bin',ha='center')

    sca(axL[2,1])
    ylabel('Cumulative Fraction of Candidates',ha='center')

    fig.subplots_adjust(wspace=0.3,top=0.99,bottom=0.07,left=0.1,right=0.9)

def prad_hist_stacked_small(prad_bins=[0.5,1.0,2.0,4.0,8.0,16]):
    lamo = cksmet.io.load_table('lamost-dr2',cache=1)
    cks =  cksphys.io.load_table(
        'cks+nea+iso-floor+huber-phot+furlan',cache=1,
        cachefn='../CKS-Physical/load_table_cache.hdf'
    )
    cks = add_cuts(cks)

    smet_bins = arange(-0.8,0.61,0.1)
    nrows = len(prad_bins) -1 
    ncols = 2
    width = 2 * ncols
    height = 2 * nrows
    fig, axL = subplots(
        nrows=nrows,ncols=ncols,figsize=(width,height),sharex=True,sharey=True
    )

    def histstack(df, axL):
        for i in range(nrows):
            sca(axL[nrows-1-i])
            prad1 = prad_bins[i]
            prad2 = prad_bins[i+1]
            cut = df[df.iso_prad.between(prad1,prad2)]
            hist(lamo.lamo_smet, bins=smet_bins, normed=True,histtype='step',lw=2)
            hist(cut.cks_smet, bins=smet_bins, normed=True,histtype='step',lw=2)
            gca().xaxis.set_major_locator(MaxNLocator(4))
            xlim(-0.8,0.6)
            
    # Lamost taps out at kepmag = 14, so let's restrict
    i = 0 
    cut = cks
    histstack(cut,axL[:,i])
    sca(axL[0,i])
    title('Full ({})'.format(len(cut)))

    for j in range(nrows):
        sca(axL[nrows-1-j,i])
        prad1 = prad_bins[j]
        prad2 = prad_bins[j+1]
        ylabel('$R_p = {:.1f}-{:.1f}$'.format(prad1,prad2))
    axL[-1,0].set_xlabel('[Fe/H]')

    i+=1

    cut = cks[cks['isfaint'.split()].sum(axis=1)==0]
    histstack(cut,axL[:,i])
    sca(axL[0,i])
    title('Kp < 14.2 ({})'.format(len(cut))) 
    i+=1

    cut = cks[cks['isfaint isneafp isfp isbadprad isdiluted isgrazing islongper issubgiant'.split()].sum(axis=1)==0]
    histstack(cut,axL[:,i])
    sca(axL[0,i])
    title('All Cuts ({})'.format(len(cut))) 
    i+=1

    fig.set_tight_layout(True)

import string

def cuts():
    df =  cksphys.io.load_table(
        'cks+nea+iso-floor+huber-phot+furlan',cache=1,
        cachefn='../CKS-Physical/load_table_cache.hdf'
    )

    cuttypes = cksmet.cuts.plnt_cuttypes
    nrows = len(cuttypes)

    nrows = 3
    ncols = 3

    width = ncols * 2.5
    height = nrows * 2.5

    fig, axL = subplots(
        nrows=nrows,ncols=ncols,figsize=(width,height),
        sharex=True,sharey=True
    )
    axL = axL.flatten()

    bpass = np.zeros(len(df))
    iplot = 0

    for cuttype in cuttypes:
        ax = axL[iplot] 
        sca(ax)
        
        key = 'is'+cuttype
        obj = cksmet.cuts.get_cut(cuttype)
        cut = obj(df,'cks')
        df[key] = cut.cut()
        bpass += df[key].astype(int)
        plotkw = dict(ms=4)
        plot(df.iso_prad, df.cks_smet,'.',color='LightGray',**plotkw)
        dfcut = df[bpass==0]
        plot(dfcut.iso_prad, dfcut.cks_smet,'.',**plotkw)

        
        _text = cut.plotstr + ' ({})'.format(len(dfcut))
        textkw = dict(fontsize='small',transform=ax.transAxes, ha='right')
        text(0.95,0.05, _text, **textkw)

        semilogx()
        xlim(0.5,32)
        xt = [0.5,1,2,4,8, 16, 32]
        xticks(xt,xt)
        iplot+=1

    axL = axL.reshape(nrows,ncols)
    setp(axL[-1,:],xlabel='$R_p\, (R_{\oplus})$')
    setp(axL[:,0],ylabel='[Fe/H]')
    fig.set_tight_layout(True)

    for ax, letter in zip(axL.flatten(),string.ascii_lowercase):
        sca(ax)
        add_anchored(
            letter,loc=2, frameon=False, 
            prop=dict(size='large', weight='bold')
        )
