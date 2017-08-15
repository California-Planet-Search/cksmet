from matplotlib.pylab import *
import seaborn as sns
import matplotlib.patheffects as path_effects
from matplotlib.transforms import blended_transform_factory as btf
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import pandas as pd 

import cksmet.analysis
import cksmet.occur

sns.set_style('ticks')
sns.set_color_codes()

def fig_checkerboard():
    occ  = cksmet.io.load_object('occur-nsmet=1') 

    prob_det_min = 0.1
    epsilon = 0.0001

    fig,ax = subplots(figsize=(8,6))

    # Completeness contour
    # levels = [prob_det_min,prob_det_min + epsilon]
    # contour(dsf.perc, dsf.pradc, dsf.prob_det, levels=levels,color='k')
    ds = occ.df.groupby(['perc','pradc']).first().to_xarray()
    ds = ds.transpose('pradc','perc')

    per = np.hstack([np.array(ds.per1[0]),np.array(ds.per2[0,-1])])
    prad = np.hstack([np.array(ds.prad1[:,0]),np.array(ds.prad2[-1,0])])
    X,Y = meshgrid(per,prad)

    ds = ds.where((ds.per2 < 350) & (ds.prob_det_mean > 0.25))
    ar = ma.masked_invalid(np.array(ds.rate))
    ar = log10(ar)
    pcolormesh(X,Y,ar,cmap='YlGn',vmin=-4,vmax=-1)
    #pcolormesh(X,Y,ar,cmap='YlGn',vmin=-4,vmax=-1)
    colorbar(ticks=[-4,-3,-2,-1])

    plot(occ.plnt.per,occ.plnt.prad,'.',ms=6,color='Tomato')
    df = ds.to_dataframe()
    annotate_checkerboard(df)
    loglog()
    label_checkerboard()
    fig.set_tight_layout(True)
    return occ

def label_checkerboard():
    fig = gcf()
    ax = gca()
    yt = [0.5, 1, 2, 4, 8, 16, 32]
    xt = [0.3,1, 3, 10, 30, 100, 300]
    xticks(xt,xt)
    yticks(yt,yt)
    xlim(0.1, 1000)
    ylim(0.3, 40)
    xlabel('Orbital Period (days)')
    ylabel('Planet Size (Earth-radii)')
    minorticks_off()
    fig.set_tight_layout(True)

def annotate_checkerboard(df,print_ul=False):
    ax = gca()
    for i, row in df.iterrows():
        xy = np.array([row.per1, row.prad1])
        width = row.per2 - row.per1
        height = row.prad2 - row.prad1
        if ~np.isnan(row.rate_ul):
            if print_ul:
                txt = "<{:.2f}".format(row.rate_ul * 100)
            else:
                txt=""
            prate = 0
        else:
            prate = row.rate * 100
            txt = "{:.2f}".format(prate)

        at =  matplotlib.patches.Rectangle(
            xy, width, height,ec='LightGray',lw=0.5, fc='none'
        )

        ax.add_artist(at)

        text = ax.annotate(txt, xy=xy, xytext=xy*1.05,size='x-small')
        text.set_path_effects(
            [
                path_effects.Stroke(linewidth=1.5, foreground='white'),
                path_effects.Normal()
            ]
        )

    '''
    for key in regions.keys():
        box = regions[key]
        xy = np.array([box[0][0],box[1][0]])
        width = box[0][1] - box[0][0]
        height = box[1][1] - box[1][0]
        at =  matplotlib.patches.Rectangle(
            xy, width, height,ec='k',lw=2, fc='none'
        )
        ax.add_artist(at)
    '''
# Occurrence period

def per(occur,**kw):
    yerr = np.array(occur['rate_err2 rate_err1'.split()]).T
    yerr[0] *= -1 
    semilogy()
    errorbar(occur.perc,occur.rate,yerr=yerr,fmt='o',**kw)
    if kw.has_key('label'):
        kw['label'] = ''
    plot(occur.perc,occur.rate_ul,marker='v',lw=0,**kw)

def per_sample(key,**kw):
    semilogx()
    fit = cksmet.io.load_object(key)
    per(fit.occur,**kw)
    peri = np.logspace(log10(1),log10(350), 100)
    fit_samples = fit.sample(peri,1000)
    fit_best = fit.model(fit.pfit,peri)
    p16, p50, p84 = np.percentile(fit_samples,[16,50,84],axis=0)
    fill_between(peri, p16, p84, color='pink',zorder=1)
    plot(peri,fit_best,color='r',alpha=1)

def fig_per_small2():
    fig,axL = subplots(ncols=2,figsize=(8.5,4),sharex=True,sharey=True)
    sca(axL[0])
    per_sample('fit_per-sub-se',color='b',label='< 0')
    per_sample('fit_per-sup-se',color='g',label='> 0')
    title('Super-Earths')
    legend(title='[Fe/H]',frameon=True)
    fig_label("a")

    sca(axL[1])
    per_sample('fit_per-sub-sn',color='b',label='< 0')
    per_sample('fit_per-sup-sn',color='g',label='> 0')
    title('Sub-Neptunes')
    legend(title='[Fe/H]',frameon=True)
    fig_label("b")

    xt = [1,3,10,30,100,300]
    xticks(xt,xt)
    setp(axL, xlabel='Period (days)', ylim=(3e-4,1), xlim = (1,300))
    setp(axL[0], ylabel='Planets Per Star')
    fig.set_tight_layout(True)

def smet(occur,**kw):
    yerr = np.array(occur['rate_err2 rate_err1'.split()]).T
    yerr[0] *= -1 
    semilogy()
    errorbar(occur.smetc,occur.rate,yerr=yerr,fmt='o',**kw)
    if kw.has_key('label'):
        kw['label'] = ''
    plot(occur.smetc,occur.rate_ul,marker='v',lw=0,**kw)

def smet_sample(key, **kw):
    fit = cksmet.io.load_object(key)
    smet(fit.occur,**kw)
    smeti = np.linspace(-0.4,0.4,100)
    fit_samples = fit.sample(smeti,1000)
    fit_best = fit.model(fit.pfit,smeti)
    p16, p50, p84 = np.percentile(fit_samples,[16,50,84],axis=0)
    fill_between(smeti, p16, p84, color='pink',zorder=1)
    plot(smeti,fit_best,color='r',alpha=1)

def fig_smet_small4():
    fig, axL = subplots(ncols=2,nrows=2,figsize=(8.5,7),sharex=True)
    sca(axL[0,0]) 
    smet_sample('fit_smet-hot-se',color='b',)
    title('Hot Super-Earths')
    fig_label("a")

    sca(axL[0,1]) 
    smet_sample('fit_smet-hot-sn',color='b',)
    title('Hot Sub-Neptunes')
    fig_label("b")

    sca(axL[1,0]) 
    smet_sample('fit_smet-warm-se',color='b',)
    title('Warm Super-Earths')
    fig_label("c")

    sca(axL[1,1]) 
    smet_sample('fit_smet-warm-sn',color='b',)
    title('Warm Sub-Neptunes')
    fig_label("d")

    setp(axL[1,:],xlabel='[Fe/H] (dex)')
    setp(axL[:,0],ylabel='Planets Per Star')
    
    setp(axL[1,:],ylim=(1e-1,6e-1))
    setp(axL[0,:],ylim=(1e-3,1e-1))
    fig.set_tight_layout(True)
    xlim(-0.4,0.4)

def fig_smet_large4():
    fig, axL = subplots(ncols=2,nrows=2,figsize=(8.5,7),sharex=True)

    df = cksmet.io.load_table('occur-nper=2-nsmet=5',cache=1)
    df = df[df.smetc.between(-0.4,0.4)]

    sca(axL[0,0]) 
    cut = df.ix['hot','ss']
    smet(cut,color='b',)
    title('Hot Sub-Saturns')
    fig_label("a")

    sca(axL[0,1]) 
    smet_sample('fit_smet-warm-ss',color='b')
    title('Warm Sub-Saturns')
    fig_label("b")

    sca(axL[1,0]) 
    smet_sample('fit_smet-hot-jup',color='b',)
    title('Hot Jupiters')
    fig_label("c")

    sca(axL[1,1]) 
    cut = df.ix['warm','jup']
    smet(cut,color='b')
    title('Warm Jupiters')
    fig_label("d")

    setp(axL[1,:],xlabel='[Fe/H] (dex)')
    setp(axL[:,0],ylabel='Planets Per Star')
    
    setp(axL[0,:],ylim=(3e-4,3e-1))
    setp(axL[1,:],ylim=(1e-4,1e-1))
    fig.set_tight_layout(True)
    xlim(-0.4,0.4)

def add_anchored(*args,**kwargs):
    """
    Parameters
    ----------
    s : string
        Text.

    loc : str
        Location code.

    pad : float, optional
        Pad between the text and the frame as fraction of the font
        size.

    borderpad : float, optional
        Pad between the frame and the axes (or *bbox_to_anchor*).

    prop : `matplotlib.font_manager.FontProperties`
        Font properties.
    """

    bbox = {}
    if kwargs.has_key('bbox'):
        bbox = kwargs.pop('bbox')
    at = AnchoredText(*args, **kwargs)
    if len(bbox.keys())>0:
        plt.setp(at.patch,**bbox)

    ax = plt.gca()
    ax.add_artist(at)

def fig_label(text):
    add_anchored(
        text, loc=2, frameon=True, 
        prop=dict(size='large', weight='bold'),
        bbox=dict(ec='none', fc='w', alpha=0.8)
    )
