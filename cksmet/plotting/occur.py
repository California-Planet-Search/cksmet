from matplotlib.pylab import *
import xarray as xr
import seaborn as sns
sns.set_style('whitegrid')

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
