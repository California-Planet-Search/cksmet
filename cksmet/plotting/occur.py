from matplotlib.pylab import *
import xarray as xr
import seaborn as sns
sns.set_style('ticks')
import matplotlib.patheffects as path_effects
import cksmet.analysis
import cksmet.occur
from cksmet.grid import *
from matplotlib.transforms import blended_transform_factory as btf

def histograms(query):
    smetbins = [
        [-0.75,-0.45],
        [-0.45,-0.15],
        [-0.15,0.15],
        [0.15,0.45]
    ]
    ydata_nstars = 0.9
    ydata_occur = 0.85
    fig,ax = subplots(figsize=(6,4))

    for i in range(len(smetbins)):
        smetbin = smetbins[i]
        cachefn = cksmet.analysis.cachefn(smetbin)
        occ  = cksmet.analysis.load_occur(
            cache=1, cachefn=cachefn, smetbin=smetbin
        ) 
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
            texoccur = "${:.1f}^{{{:+.1f}}}_{{{:+.1f}}}\%$".format(
                prate, prate_err1, prate_err2
            )
            
            y = [prate]
            yerr = [[-prate_err2],[prate_err1]]
            errorbar([x],[y], xerr=xerr, yerr=yerr,**errkw)
        else:
            out = cksmet.occur.combine_cells(cut)
            prate_ul = out['rate_ul'] * 100
            soccur = "< {:.2f}%".format(prate_ul)
            texoccur = "$< {:.1f}\%$".format(prate_ul)

            errorbar(
                [x],[prate_ul], xerr=xerr, yerr=[prate_ul], uplims=True, 
                capsize=5,**errkw
            )
        
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



regions = {
    'HotJ':[[1,10],[8,22.6]],
    'HotSupE':[[3,30],[1,2]],
    'HotSubN':[[3,30],[2,4]],
    'CoolSubN':[[30,300],[2,4]],
}


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

def mesfac(st):
    fig,axL = subplots(ncols=2,figsize=(8,4))
    sca(axL[0])
    loglog()
    plot(st.plnt_mes.mes,st.plnt_mes.mes_formula,'.')
    xlabel('MES (TPS)')
    ylabel('MES (Formula)')
    s = '{:.2f} $\pm$ {:.2f} (dex)'.format(
        st.logmes_diff_std, st.logmes_diff_med
    )
    text(0.6,0.1,s,transform=axL[0].transAxes)
    sca(axL[1])
    loglog()
    plot(st.plnt_mes.mes,st.plnt_mes.mes_formula_scaled,'.')
    xlabel('MES (TPS)')
    ylabel('MES (Formula-Scaled)')
    x = [10,1e4]
    for ax in axL:   
        sca(ax)
        ax.plot(x,x)
        xlim(*x)
        ylim(*x)
        
    fig.set_tight_layout(True)

def compare_prob_det_direct_interp(comp):
    df = comp.grid.ds.to_dataframe()
    df['prob_det_direct'] = vstack([comp.prob_det(row.per1,row.prad1,method='direct') for i,row in df.iterrows()])
    df['prob_det_interp'] = vstack([comp.prob_det(row.per1,row.prad1,method='interp') for i,row in df.iterrows()])
    ds = df.to_xarray()
    fig,axL = subplots(ncols=3,figsize=(12,4))
    sca(axL[0])
    loglog()
    ds.prob_det_direct.plot.contourf(y='prad',x='per',ax=axL[0],fig=fig)

    sca(axL[1])
    loglog()
    ds.prob_det_interp.plot.contourf(y='prad',x='per',ax=axL[1],fig=fig)
    
    sca(axL[2])
    loglog()
    rat = ds.prob_det_direct / ds.prob_det_interp
    rat.plot.contourf(y='prad',x='per',ax=axL[2],fig=fig, vmax=1.01,vmin=0.99)
    fig.set_tight_layout(True)

def label():
    xlabel('Orbital Period (days)')
    ylabel('Planet Size (Earth-radii)')
    xt = [1,3,10,30,100,300]
    yt = [0.5,1,2,4,8,16]
    xticks(xt,xt)
    yticks(yt,yt)
    fig = gcf()
    rename()

def rename():
    """
    Plotting routines will give short machine readable names, here we
    change to human-readable names

    """
    fig = gcf()
    namemap = {
        'prob_det':'Prob Det.',
        'prob_tr':'Prob Tr.',
        'prob_trdet':'Prob Tr. and Det.',
        'plnt_occur':'Planet Occurrence per Bin',
        'prad':'Planet Size (Earth-radii)',
        'fstars':'Completeness',
    }
    for o in fig.findobj(matplotlib.text.Text):
        for oldname, newname in namemap.iteritems():
            if o.get_text()==oldname:
                o.set_text(newname)

def plot_completeness(comp):
    ds = comp.grid.ds
    ds['prob_det'] = comp.grid.ds['prob_det']
    ds['prob_trdet'] = comp.grid.ds['prob_tr'] * ds['prob_det']

    kw = dict(
        x='per', y='prad',cmap=cm.Blues_r,vmax=1,levels=arange(0,1.1,0.1),
        zorder=0
    )

    fig,axL = subplots(ncols=2,figsize=(10,4))
    sca(axL[0])
    loglog()
    ds['prob_det'].plot.contourf(**kw)
    label()

    sca(axL[1])
    loglog()
    kw['levels'] = [0.000,0.001,0.010,0.020,0.050,0.100]
    ds['prob_trdet'].plot.contourf(**kw)
    fig = gcf()
    fig.set_tight_layout(True)
    label()
    ylabel('')

smetbins = [
    [-0.70, -0.45],
    [-0.45, -0.15],
    [-0.15, 0.15],
    [0.15,  0.45],
]


def occur_smet_fourpanel():
    """

    """

    test = xr.open_dataset('data/test.nc') 
    fig,axL = subplots(nrows=2,ncols=2)
    axL = axL.flatten()
    
    for i, smetbin in enumerate(smetbins):
        sca(axL[i])
        loglog()
        
        smet0 = smetbin[0]
        smet1 = smetbin[1]
        dscut = test.where(( smet0 < test.smet) & (test.smet < smet1))
        dscut2d = dscut.sum(dim=['smet'])
        dscut2d.plnt_occur.plot.imshow(vmax=0.03)

        _title = '[Fe/H] = $\mathregular{{ {0:+.2f} - {1:+.2f} }}$'.format(*smetbin)
        title(_title)
    
    setp(axL,xlim=(1,300))

    yt = [0.5, 1, 2, 4, 8, 16, 32]
    xt = [1, 3, 10, 30, 100, 300]
    for ax in axL:
        sca(ax)
        yticks(yt,yt)
        xticks(xt,xt)
        
    fig.set_tight_layout(True)
