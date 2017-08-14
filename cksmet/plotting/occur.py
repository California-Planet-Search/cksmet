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

def histograms(smetbins, occL, query, tworuns=True):
    if tworuns==True:
        minrate = histograms(smetbins, occL, query,tworuns=False),
        yl = list(ylim())
        yuplim = 0.1 * (yl[1] - yl[0])
        cla() 
    else:
        yuplim=nan
        fig, ax = subplots(figsize=(6,4))

    ax = gca()
    ydata_nstars = 0.9
    ydata_occur = 0.85
    minrate = 1
    for i in range(len(smetbins)):
        smetbin = smetbins[i]
        occ = occL[i]

        cut = occ.grid_loguni.ds.to_dataframe().query(query)
        x = 0.5 * (smetbin[0] + smetbin[1])
        xerr = 0.5 * (smetbin[1] - smetbin[0])
        errkw = dict(color='k')

        sferng = "[Fe/H] = [{:+.2f},{:+.2f}]".format(*smetbin)

        trans = btf(ax.transData, ax.transAxes)
        nplanets = cut.nplanets.sum()
        if nplanets > 0:
            out = cksmet.occur.combine_cells(cut)
            prate = out['rate'] * 100
            prate_err1 = out['rate_err1'] * 100
            prate_err2 = out['rate_err2'] * 100

            soccur = "{:.2f} {:+.2f}/{:+.2f}%".format(
                prate, prate_err1, prate_err2
            )
            texoccur = "${:.2f}^{{{:+.2f}}}_{{{:+.2f}}}\%$".format(
                prate, prate_err1, prate_err2
            )
            
            y = [prate]
            yerr = [[-prate_err2],[prate_err1]]
            errorbar([x],[y], xerr=xerr, yerr=yerr,**errkw)
            minrate = min(minrate,prate)
        else:
            out = cksmet.occur.combine_cells(cut)
            prate_ul = out['rate_ul'] * 100
            soccur = "< {:.2f}%".format(prate_ul)
            texoccur = "$< {:.2f}\%$".format(prate_ul)

            yuplim = min(yuplim,prate_ul)
            errorbar(
                [x],[prate_ul], xerr=xerr, yerr=[yuplim], uplims=True, 
                capsize=5,**errkw
            )
            minrate = min(minrate,prate_ul)
        
        print "{}, {}".format(sferng, soccur)
        print "Mean ntrials {:.1f}".format(out['ntrials'])
        print "" 
        text(smetbin[0], ydata_occur, texoccur, transform=trans)
        
        texnstars = "${:.0f}/{:.0f}$".format(nplanets, occ.nstars)
        text(smetbin[0], ydata_nstars, texnstars, transform=trans)

    xdata_label = smetbins[0][0] - 0.2
    text(xdata_label, ydata_nstars, " $N_p/N_\star$", transform=trans)
    text(xdata_label, ydata_occur, " $f$", transform=trans)

    yl = list(ylim())
    yl[0] = 0
    yl[1] *= 1.2
    ylim(*yl)

    xl = list(xlim())
    xl[0] = xdata_label
    xlim(*xl)
    xlabel('[Fe/H]')
    ylabel('Planets per 100 Stars')
    return minrate 

def checkerboard(smetbin, print_ul=False):
    cachefn = cksmet.analysis.cachefn(smetbin)
    occ  = cksmet.analysis.load_occur(cache=1,cachefn=cachefn,smetbin=smetbin) 
    #downsamp = {'per':2,'prad':4}
    #occ.set_grid_loguni(downsamp)
    prob_det_min = 0.1
    epsilon = 0.0001

    fig,ax = subplots(figsize=(8,6))
    plot(occ.plnt.per,occ.plnt.prad,'.',ms=6,color='Tomato')
    dsf = occ.grid_fine.ds
    
    levels = [prob_det_min,prob_det_min + epsilon]
    contour(dsf.perc, dsf.pradc, dsf.prob_det, levels=levels,color='k')

    ds = occ.grid_loguni.ds.transpose('prad','per')
    per = np.hstack([np.array(ds.per1[0]),np.array(ds.per2[0,-1])])
    prad = np.hstack([np.array(ds.prad1[:,0]),np.array(ds.prad2[-1,0])])
    X,Y = meshgrid(per,prad)

    ds = ds.where(
        (ds.per1 > 0.3) & (ds.per2 < 350) & (ds.prob_det_min > prob_det_min) 
    )
    ar = ma.masked_invalid(np.array(ds.rate))
    ar = log10(ar)
    pcolormesh(X,Y,ar,cmap='YlGn',vmin=-4,vmax=-1)
    colorbar(ticks=[-4,-3,-2,-1])

    df = ds.to_dataframe()
    annotate_checkerboard(df)
    loglog()
    label_checkerboard()
    return occ

def checkerboard_ratio(smetbin, print_ul=False):
    cachefn = cksmet.analysis.cachefn([-0.75,0.5])
    occ0  = cksmet.analysis.load_occur(cache=1,cachefn=cachefn,smetbin=smetbin)
    cachefn = cksmet.analysis.cachefn(smetbin)
    occ  = cksmet.analysis.load_occur(cache=1,cachefn=cachefn,smetbin=smetbin) 

    downsamp = {'per':2,'prad':4}
    occ0.set_grid_loguni(downsamp)
    occ.set_grid_loguni(downsamp)
    prob_det_min = 0.25 

    fig,ax = subplots(figsize=(8,6))
    loglog()
    plot(occ.plnt.per,occ.plnt.prad,'.',ms=6,color='Tomato')
    dsf = occ.grid_fine.ds
    contour(dsf.perc,dsf.pradc,dsf.prob_det,levels=[0.25,0.251],color='k')

    ds0 = occ0.grid_loguni.ds.transpose('prad','per')
    ds = occ.grid_loguni.ds.transpose('prad','per')
    per = np.hstack([np.array(ds.per1[0]),np.array(ds.per2[0,-1])])
    prad = np.hstack([np.array(ds.prad1[:,0]),np.array(ds.prad2[-1,0])])
    X,Y = meshgrid(per,prad)

    ds = ds.where(
        (ds.per1 > 0.3) & (ds.per2 < 350) & (ds.prob_det_min > 0.25) 
    )
    #ratio
    ds['lograte_diff'] = np.log10(ds.rate) - np.log10(ds0.rate)
    ar = ma.masked_invalid(ds.lograte_diff)

    #ar = log10(ar)
    pcolormesh(X,Y,ar,cmap='coolwarm',vmin=-1,vmax=1)
    colorbar(ticks=[-2,-1,0,1,2])

    df = ds.to_dataframe()
    annotate_checkerboard_diff(df)
    label_checkerboard()
    return ds0,ds

def checkerboard_ratio(smetbin, print_ul=False):
    cachefn = cksmet.analysis.cachefn([-0.75,0.5])
    occ0  = cksmet.analysis.load_occur(cache=1,cachefn=cachefn,smetbin=smetbin)
    cachefn = cksmet.analysis.cachefn(smetbin)
    occ  = cksmet.analysis.load_occur(cache=1,cachefn=cachefn,smetbin=smetbin) 

    downsamp = {'per':2,'prad':4}
    occ0.set_grid_loguni(downsamp)
    occ.set_grid_loguni(downsamp)
    prob_det_min = 0.25 

    fig,ax = subplots(figsize=(8,6))
    loglog()
    plot(occ.plnt.per,occ.plnt.prad,'.',ms=6,color='Tomato')
    dsf = occ.grid_fine.ds
    contour(dsf.perc,dsf.pradc,dsf.prob_det,levels=[0.25,0.251],color='k')

    ds0 = occ0.grid_loguni.ds.transpose('prad','per')
    ds = occ.grid_loguni.ds.transpose('prad','per')
    per = np.hstack([np.array(ds.per1[0]),np.array(ds.per2[0,-1])])
    prad = np.hstack([np.array(ds.prad1[:,0]),np.array(ds.prad2[-1,0])])
    X,Y = meshgrid(per,prad)

    ds = ds.where(
        (ds.per1 > 0.3) & (ds.per2 < 350) & (ds.prob_det_min > 0.25) 
    )
    #ratio
    ds['lograte_diff'] = np.log10(ds.rate) - np.log10(ds0.rate)
    ar = ma.masked_invalid(ds.lograte_diff)

    #ar = log10(ar)
    pcolormesh(X,Y,ar,cmap='coolwarm',vmin=-1,vmax=1)
    colorbar(ticks=[-2,-1,0,1,2])

    df = ds.to_dataframe()
    annotate_checkerboard_diff(df)
    label_checkerboard()
    return ds0,ds

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

    for key in regions.keys():
        box = regions[key]
        xy = np.array([box[0][0],box[1][0]])
        width = box[0][1] - box[0][0]
        height = box[1][1] - box[1][0]
        at =  matplotlib.patches.Rectangle(
            xy, width, height,ec='k',lw=2, fc='none'
        )
        ax.add_artist(at)


def annotate_checkerboard_diff(df,print_ul=False):
    ax = gca()
    for i, row in df.iterrows():
        xy = np.array([row.per1, row.prad1])
        width = row.per2 - row.per1
        height = row.prad2 - row.prad1
        if ~np.isnan(row.rate_ul):
            txt=""
        else:
            prate = row.lograte_diff 
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



def binned_period(df, **kwargs):
    kw = dict(fmt='_',mew=2,ms=7,capsize=2,**kwargs)
    i_count = 0 
    for i, row in df.iterrows():
        if i_count!=0:
            kw['label']=None
        i_count+=1

        rate, rate_err1, rate_err2, rate_ul, stats = \
            cksmet.occur.binomial_rate(row.ntrial, row.nplnt)
        #print row.nplnt, row.ntrial,rate, rate_err1, rate_err2, rate_ul
        if row.nplnt > 1:
            yerr = [[-rate_err1],[rate_err2]] 
            errorbar([row.perc],rate,yerr=yerr,**kw)
        else:
            errorbar([row.perc],rate_ul,yerr=rate_ul*0.2,uplims=True,**kw)
            
    xlabel('Orbital Period (days)')
    ylabel('Planets Per Star')


'''
axvline(res.params['per0'],color=color,ls='--',lw=1)
ax = gca()
trans = btf(ax.transData, ax.transAxes) 
s = "P = {:.1f} days".format(res.params['per0'].value)
text(res.params['per0'],0.925,s, color=color,rotation=90, transform=trans, ha=ha,size='medium')
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
    fit = cksmet.io.load_fit(key)
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
    fit = cksmet.io.load_fit(key)
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
