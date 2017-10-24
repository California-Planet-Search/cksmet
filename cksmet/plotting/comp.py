from matplotlib.pylab import *
import seaborn as sns
from cksmet.plotting.config import *
import cksmet.io

def label():
    xlabel('Orbital Period (days)')
    ylabel('Planet Size (Earth-radii)')
    xt = [1,3,10,30,100,300]
    yt = [0.5,1,2,4,8,16]
    xticks(xt,xt)
    yticks(yt,yt)
    fig = gcf()
    rename()

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

def fig_prob_detect_transit():
    sns.set_style('whitegrid')
    sns.set_context('paper')
    fig,axL = subplots(ncols=2,figsize=(7.25,2.75),)

    comp = cksmet.io.load_object('comp',cache=1)

    ds = comp.grid.ds
    ds['prob_det'] = comp.grid.ds['prob_det']
    ds['prob_trdet'] = comp.grid.ds['prob_tr'] * ds['prob_det']

    kw = dict(
        x='per', y='prad',cmap=cm.Blues_r,vmax=1,levels=arange(0,1.1,0.1),
        zorder=0,add_colorbar=False
    )

    sca(axL[0])
    loglog()
    im = ds['prob_det'].plot.contourf(**kw)
    cbar = colorbar(im,shrink=0.5)
    cbar.set_label(r'$\left< p_{\mathrm{det}} \right>$',size='small')
    t = cbar.ax.get_yticklabels()
    setp(t,size='x-small')

    label()
    sca(axL[1])
    loglog()
    kw['levels'] = [0.000,0.001,0.010,0.020,0.050,0.100]
    im = ds['prob_trdet'].plot.contourf(**kw)
    cbar = colorbar(im,shrink=0.5)
    cbar.set_label(r'$\left< p_{\mathrm{tr}} p_{\mathrm{det}} \right>$',size='small')
    t = cbar.ax.get_yticklabels()
    setp(t,size='x-small')
    label()

    setp(axL,xlim=(1,300),ylim=(0.5,32))
    tight_layout(rect=[0.05,0.01,0.99,0.99],pad=0.05)
