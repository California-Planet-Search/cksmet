from matplotlib.pylab import *
import seaborn as sns
import matplotlib.patheffects as path_effects
from matplotlib.transforms import blended_transform_factory as btf
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import pandas as pd 
import cksmet.io
from scipy import ndimage as nd

import cksmet.analysis
import cksmet.occur
from matplotlib.patches import Rectangle
sns.set_style('ticks')
#sns.set_palette('bright')
sns.set_color_codes()

ptcolor = {
    'se':'g',
    'sn':'b',
    'ss':'y',
    'jup':'r'
}

bdcolor = {
    'se':'light green',
    'sn':'light blue',
    'ss':'light mustard',
    'jup':'light pink'
}

for k in ptcolor.keys():
    bdcolor[k] = sns.xkcd_rgb[bdcolor[k]]

namedict = {
    'se':'Super-Earths',
    'sn':'Sub-Neptunes',
    'ss':'Sub-Saturns',
    'jup':'Jupiters'
}

sizedict = {
    'se':'$R_P$ = 1.0$-$1.7 $R_E$',
    'sn':'$R_P$ = 1.7$-$4.0 $R_E$',
    'ss':'$R_P$ = 4.0$-$8.0 $R_E$',
    'jup':'$R_P$ = 8.0$-$24.0 $R_E$'
}

bin_per = 0.25 # period bin width in dex
bin_prad = 0.1 # prad bin width in dex

def calc_contour():
    per1 = 0.5
    per2 = 1000
    prad1 = 0.25
    prad2 = 32
    ssamp_per = 10 # supersample factor
    ssamp_prad = 10 # supersample factor
    eps = 1e-3
    smet_bins = [-1,0.5]
    df = []
    count = 0 
    for i in range(ssamp_per):
        for j in range(ssamp_prad):
            if count < 2:
                verbose=1
            else:
                verbose=0
            shift_logper = bin_per / ssamp_per * i
            shift_logprad = bin_prad / ssamp_prad * j
            logper1 = np.log10(per1) + shift_logper
            logper2 = np.log10(per2) + shift_logper
            logprad1 = np.log10(prad1) + shift_logprad
            logprad2 = np.log10(prad2) + shift_logprad
            per_bins = 10**(np.arange(logper1,logper2+eps,bin_per))
            prad_bins = 10**(np.arange(logprad1,logprad2+eps,bin_prad))
            print per_bins
            print prad_bins

            occ = cksmet.analysis.compute_binned_occurrence(
                per_bins, prad_bins, smet_bins, verbose=verbose
            )
            occ.df['count'] = count
            df.append(occ.df)
            count+=1
            
    df = pd.concat(df)
    df.to_hdf('test.hdf','test')
    return df

def fig_contour_linear(scale='linear'):
    contour(scale='linear')

def fig_contour_log(scale='log'):
    contour(scale='log')

def fig_contour_all(scale='log'):
    import seaborn as sns
    sns.set_context('paper',font_scale=1.0)
    sns.set_style('ticks')
    fig, axL = subplots(nrows=2,ncols=2,figsize=(8,6.5),sharex=False,sharey=False)

    sca(axL[0,0])
    fig_label("a")
    contour(scale='linear',draw_colorbar=False,plot_interval=True)
    
    sca(axL[0,1])
    fig_label("b")
    caxheight = 0.25
    caxcenter = 0.8
    caxleft = 0.92
    cax = plt.axes([caxleft, caxcenter - 0.5*caxheight, .01, caxheight])
    sca(axL[0,1])    
    contour(scale='linear',plot_planetpoints=False, draw_colorbar=True,cax=cax,label=True)

    sca(axL[1,0])
    fig_label("c")
    contour(scale='log',draw_colorbar=False)
    
    sca(axL[1,1])
    fig_label("d")
    caxcenter = 0.3
    cax = plt.axes([caxleft, caxcenter - 0.5*caxheight, .01, caxheight])
    sca(axL[1,1])    
    contour(scale='log',plot_planetpoints=False, draw_colorbar=True,cax=cax,label=True)

    tight_layout(rect=[0.01,0.01,caxleft-0.01,0.99],pad=0.00)
    setp(axL[0,:],xlabel='')
    setp(axL[:,1],ylabel='')
    

def contour(scale='linear', plot_planetpoints=True, plot_interval=False, 
            draw_colorbar=True,cax=None,plot_completeness=True,label=False):

    ax = gca()
    tax = gca().transAxes

    binsize_prad = 0.1 # dex
    binsize_per = 0.5 # dex

    test = pd.read_hdf('test.hdf','test')
    temp = test.query('per2 < 350 and prob_det_mean > 0.25').dropna(subset=['rate']).sort_values(by='rate')
    testmin = temp.iloc[0]
    testmax = temp.iloc[-1]
    print "max {perc} {pradc} {rate}".format(**testmax)
    print "min {perc} {pradc} {rate}".format(**testmin)
    test = test.groupby(['per1','prad1']).first()

    ds = test.to_xarray()
    cut = ds
    #cut = ds.where((ds.per2 < 350) & (ds.prob_det_mean > 0.25))
    

    rate = cut.rate
    rate = rate.fillna(4e-6)
    rate = nd.gaussian_filter(rate,(2*2.5,1*2.5))

    eps = 1e-10

    X, Y = np.log10(ds.perc), np.log10(ds.pradc)

    cmap = None
    cmap = 'YlGn' #,None #'hot_r'
    #cmap = sns.light_palette((260, 75, 60), input="husl",as_cmap=True)
    #cmap = sns.cubehelix_palette(rot=-.4,as_cmap=True)
    #cmap = "Reds"
    levels = None
    cbarlabel=''
    if scale=='linear':
        levels = arange(0,5.0e-2+eps,0.005) 
        Z = rate
        kw = dict(extend='neither',cmap=cmap)
        cbarticks = levels[::2]
        cbarticklabels = ["{:.0f}".format(1e2*_yt) for _yt in cbarticks]
        kw = dict(levels=levels,extend='neither',cmap=cmap)

    # log scale to show Hot-J
    if scale=='log':
        Z = np.log10(rate)
        #levels = np.arange(-3.75,-1.5+eps,0.125) 
        #cbarticks = levels[::2]
        #cbarticklabels = ["{:.2f}".format(1e2*10**_yt) for _yt in cbarticks]
        levels = np.arange(-4,-1+eps,0.25) 
        cbarticklabels = [0.01, 0.03, 0.1, 0.3, 1, 3,10]
        cbarticks = np.log10(np.array(cbarticklabels) * 1e-2)
        kw = dict(extend='min',cmap=cmap,levels=levels,vmin=-3.99)


    cbarlabel = r"""Planets per 100 Stars per $P-R_P$ interval"""

    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    qcs = contourf(X,Y,Z, **kw)

    # plot straight contours
    #kw.pop('cmap')
    #kw.pop('extend')
    #qcs = contour(X,Y,Z,)
    #plt.clabel(qcs, inline=1, fmt='%.3f', colors='w', fontsize=1)
    if draw_colorbar:
        cbar = colorbar(qcs,cax=cax,ticks=cbarticks,)
        t = cbar.ax.set_yticklabels(cbarticklabels) 
        setp(t,size='x-small')
        cbar.set_label(cbarlabel,size='small')

    if plot_planetpoints:
        cks = cksmet.io.load_table('cks-cuts',cache=1)
        cks = cks[~cks.isany]
        x,y = log10(cks.koi_period),log10(cks.iso_prad)
        plot(x,y,'.',mfc='none',mec='red',mew=0.3,ms=5)

    # Completeness
    if plot_completeness:
        Z = np.array(ds.prob_det_mean)
        cmap = sns.light_palette("gray",as_cmap=True)
        contourf(X,Y,Z,[0,0.25],zorder=2.5,cmap=cmap,vmax=1)
        text(
            0.95,0.15,'Low Completeness',rotation=12,zorder=5,size='small',
            transform=ax.transAxes,ha='right'
        )

    if plot_interval:
        #inv = ax.transAxes.inverted()
        xyaxes = (0.6,0.9)
        xy = ax.transLimits.inverted().transform(xyaxes)
        #xy = ax.transAxes.transform((0.9,0.9))
        w = bin_per 
        h = bin_prad
        rect = Rectangle(xy, w, h,lw=1,ec='r',fc='none',zorder=4)
        ax.add_patch(rect)
        s = """\
$P-R_P$ Interval
$\Delta \log P$ = {:.2f} dex
$\Delta \log R_P$ = {:.2f} dex
""".format(w,h)
        kw = dict(
            size='x-small',zorder=5,va='top',ha='left',transform=ax.transAxes,
        )
        text(xyaxes[0]+0.07,xyaxes[1],s,**kw)


    if label:
        if scale=='linear':
            kw = dict(
                size='x-small',zorder=5,va='center',ha='center',color='red'
            )
            text(log10(30),log10(3),'Warm Sub-Neptunes', **kw)
            text(log10(20),log10(1),'Warm Super-Earths',**kw)
            text(log10(30),log10(1.7),'Radius Gap',rotation=-10,**kw)

        if scale=='log':
            kw = dict(
                size='x-small',zorder=5,va='center',ha='center',color='red'
            )
            text(log10(3),log10(20),'Hot Jupiters',**kw)
            x = [1,3.0,15,1]
            y = [1.7,1.7,10.0,10.0]
            x = np.log10(np.array(x))
            y = np.log10(np.array(y))
            plot(x,y,linestyle='--',color='red',lw=1)
            kw['va'] = 'top'
            text(log10(3),log10(10)-0.03,'Hot Planet Desert',**kw)

    xt = [1,3,10,30,100,300]
    yt = [0.5,1,2,4,8,16,32]
    xticks([log10(_xt) for _xt in xt],xt)
    yticks([log10(_yt) for _yt in yt],yt)
    xlim(log10(1),log10(300))
    ylim(log10(0.5),log10(32))

    xlabel('Orbital Period (days)')
    ylabel('Planet Size (Earth-radii)')
    
def fig_checkerboard(plot_planetpoints=True):
    occ  = cksmet.io.load_object('occur-nsmet=1') 
    sns.set_context('paper',font_scale=1.0)
    sns.set_style('ticks')
    fig,ax = subplots(figsize=(7,5.5))
    loglog()

    logper1 = np.log10(1)
    logper2 = np.log10(350)
    logprad1 = np.log10(0.5)
    logprad2 = np.log10(32)
    eps = 1e-10
    bin_per = 0.25 # period bin width in dex
    bin_prad = 0.5 * log10(2) # prad bin width in dex
    per_bins = 10**(np.arange(logper1,logper2+eps,bin_per))
    prad_bins = 10**(np.arange(logprad1,logprad2+eps,bin_prad))
    eps = 1e-3
    smet_bins = [-1,0.5]

    occ = cksmet.analysis.compute_binned_occurrence(
        per_bins, prad_bins, smet_bins
    )

    prob_det_min = 0.1
    cmap='YlGn'
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
    qcs = pcolormesh(X,Y,ar,cmap=cmap,vmin=-4,vmax=-1)

    if plot_planetpoints:
        cks = cksmet.io.load_table('cks-cuts',cache=1)
        cks = cks[~cks.isany]
        x,y = cks.koi_period,cks.iso_prad
        plot(x,y,'.',mfc='none',mec='red',mew=0.5,ms=5)

    if 1:
        caxheight = 0.5
        caxcenter = 0.5
        caxleft = 0.88
        cax = plt.axes([caxleft, caxcenter - 0.5*caxheight, .01, caxheight])
        sca(ax)    

        cbarticklabels = [0.01,0.03, 0.1, 0.3, 1, 3, 10]
        cbarticks = np.log10(np.array(cbarticklabels) * 1e-2)
        cbar = colorbar(qcs,cax=cax,ticks=cbarticks)
        #cbar = colorbar(qcs,cax=cax)
        t = cbar.ax.set_yticklabels(cbarticklabels) 
        setp(t,size='x-small')
        cbarlabel = r"""Planets per 100 Stars per Bin"""
        cbar.set_label(cbarlabel,size='small')

        #plot(occ.plnt.per,occ.plnt.prad,'o',ms=5,color='Tomato')
    df = ds.to_dataframe()
    annotate_checkerboard(df)
    label_checkerboard()
    tight_layout(rect=[0.01,0.01,0.85,0.99],pad=0)
    return occ

def label_checkerboard():
    yt = [0.5, 1, 2, 4, 8, 16, 32]
    xt = [0.3,1, 3, 10, 30, 100, 300]
    xticks(xt,xt)
    yticks(yt,yt)
    xlim(1,300)
    ylim(0.5, 32)
    xlabel('Orbital Period (days)')
    ylabel('Planet Size (Earth-radii)')
    minorticks_off()

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

def per(occur,**kw):
    yerr = np.array(occur['rate_err2 rate_err1'.split()]).T
    yerr[0] *= -1 
    semilogy()
    errorbar(occur.perc,occur.rate,yerr=yerr,fmt='o',**kw)
    if kw.has_key('label'):
        kw['label'] = ''
    plot(occur.perc,occur.rate_ul,marker='v',lw=0,**kw)

def fig_per_small2():
    fig,axL = subplots(ncols=2,figsize=(8.5,4),sharex=True,sharey=True)
    sca(axL[0])
    per_sample('fit_per-sub-se','se',color='b',label='[Fe/H] < 0')
    per_sample('fit_per-sup-se','se',color='g',label='[Fe/H] > 0')
    title('Super-Earths')
    legend(frameon=True)
    fig_label("a")

    ylabel('Planets Per 100 Stars')
    yticks_planets_per_100_stars()

    sca(axL[1])
    per_sample('fit_per-sub-sn','sn',color='b',label='[Fe/H] < 0')
    per_sample('fit_per-sup-sn','sn',color='g',label='[Fe/H] > 0')
    title('Sub-Neptunes')
    legend(frameon=True)
    fig_label("b")
    yticks_planets_per_100_stars()

    xt = [1,3,10,30,100,300]
    xticks(xt,xt)
    yticks_planets_per_100_stars()
    setp(axL, xlabel='Period (days)', ylim=(3e-4,1), xlim = (1,300))
    fig.set_tight_layout(True)


def fig_per():
    fig,axL = subplots(ncols=2,figsize=(8.5,4),sharex=True,sharey=True)
    sca(axL[0])
    per_sample('fit_per-sub-se',color='b',label='[Fe/H] < 0')
    per_sample('fit_per-sup-se',color='g',label='[Fe/H] > 0')
    title('Super-Earths')
    legend(frameon=True)
    fig_label("a")

    ylabel('Planets Per 100 Stars')
    yticks_planets_per_100_stars()

    sca(axL[1])
    per_sample('fit_per-sub-sn',color='b',label='[Fe/H] < 0')
    per_sample('fit_per-sup-sn',color='g',label='[Fe/H] > 0')
    title('Sub-Neptunes')
    legend(frameon=True)
    fig_label("b")
    yticks_planets_per_100_stars()

    xt = [1,3,10,30,100,300]
    xticks(xt,xt)
    yticks_planets_per_100_stars()
    setp(axL, xlabel='Period (days)', ylim=(3e-4,1), xlim = (1,300))
    fig.set_tight_layout(True)

def plot_rates(xk, occur, fmtkey, **kw):
    """
    Args
        xk (str): x value
        occur (pd.DataFrame): must contain rate, rate_err1, rate_err2
    """
    if xk=='smetc':
        _ptcolor = ptcolor[fmtkey]
        
    ebkw1 = dict(fmt='o',ms=6,mew=2,mfc='w',capsize=6,capthick=2,zorder=4,color=_ptcolor,**kw)
    ebkw2 = dict(**ebkw1)
    ebkw2['mfc'] = 'none'
    ebkw2['zorder'] = 5
    ebkw2['lw'] = 0
    ebkw2['zorder'] = 5

    ulkw = dict()
    ulkw['color'] = _ptcolor
    ulkw['marker'] = 'v'
    ulkw['mfc'] = 'w'
    ulkw['mew'] = 2
    ulkw['lw'] = 0
    ulkw['zorder'] = 6
    
    yerr = np.array(occur['rate_err2 rate_err1'.split()]).T
    yerr[0] *= -1 

    semilogy()

    x = occur[xk]
    y = occur.rate
    errorbar(x,y,yerr=yerr, **ebkw1)
    errorbar(x,y,yerr=yerr, **ebkw2)

    occurul = occur.dropna(subset=['rate_ul'])
    if len(occurul) >0:
        plot(x,occur.rate_ul,**ulkw)

class Sampler(object):
    def __init__(self, key, fmtkey):
        self.fit = cksmet.io.load_object(key)
        self.fmtkey = fmtkey

    def plot_band(self):
        p16, p50, p84 = np.percentile(self.fit_samples,[16,50,84],axis=0)
        _bdcolor = ptcolor[self.fmtkey]
        fill_between(self.x, p16, p84, color=_bdcolor,alpha=0.3)
        
    def plot_best(self):
        plot(self.x,self.fit_best,color=ptcolor[self.fmtkey])

    def plot_all(self):
        self.compute_samples()
        self.compute_best()
        self.plot_band()
        self.plot_best()

    def dict_to_params(self, params):
        _params = lmfit.Parameters()
        for k in params.keys():
            _params.add(k, value=params[k] )
        return _params 

    def compute_samples(self):
        psamples = self.fit.sample_chain(1000)
        fit_samples = []
        for i, params in psamples.iterrows():
            params = self.dict_to_params(params)
            fit_samples.append(self.fit.model(params,self.x))
        self.fit_samples = np.vstack(fit_samples)

    def compute_best(self):
        params = self.fit.pfit
        self.fit_best = self.fit.model(params,self.x)

import lmfit 
class SamplerSmet(Sampler):
    x = np.linspace(-0.4,0.4,10) 
    
class SamplerPerSmet(Sampler):
    smet1 = -0.4
    smet2 = 0.4
    x = linspace(smet1,smet2,10)
    smetwidth = 0.2

    def compute_samples(self):
        psamples = self.fit.sample_chain(1000)
        fit_samples = []
        for i, params in psamples.iterrows():
            fit_samples.append(self.model(params))
        self.fit_samples = np.vstack(fit_samples)

    def compute_best(self):
        p = self.fit.pfit.valuesdict()
        self.fit_best = self.model(p)

    # Integrate over the period axis
    def model(self, params):
        fit_sample = []
        for _smeti in self.x:
            _smet1 = _smeti - 0.5*self.smetwidth
            _smet2 = _smeti + 0.5*self.smetwidth
            _per1 = 1
            _per2 = 10
            lims = [[_per1,_per2],[_smet1,_smet2]]
            val = cksmet.fit.per_powerlaw_smet_exp_integral(params,lims)
            val = val / self.smetwidth # divide by the size of the bin
            fit_sample.append(val)
        fit_sample = np.array(fit_sample)
        return fit_sample


def per_sample(key, fmtkey, **kw):
    semilogx()
    fit = cksmet.io.load_object(key)
    per(fit.occur,**kw)
    peri = np.logspace(log10(1),log10(350), 100)
    fit_samples = fit.sample(peri,1000)
    fit_best = fit.model(fit.pfit,peri)
    p16, p50, p84 = np.percentile(fit_samples,[16,50,84],axis=0)
    #_bdcolor = bdcolor[fmtkey]
    #fill_between(peri, p16, p84, color=_bdcolor,zorder=1)
    fill_between(peri, p16, p84, color='pink',zorder=1)
    plot(peri,fit_best,color='r',alpha=1)


def persmet_sample(key, fmtkey, **kw):
    fit = cksmet.io.load_object(key)
    #smet(fit.occur,**kw)


def sum_cells_per(df0):
    df2 = []
    for smetc in df0.smetc.drop_duplicates():
        df = df0[df0.smetc==smetc] 
        rate = cksmet.stats.sum_cells(df.ntrial,df.nplnt)
        rate['smetc'] = smetc
        df2.append(rate)
    
    df2 = pd.DataFrame(df2)
    return df2

    
def fig_smet_warm():
    occ = cksmet.io.load_object('occur-nper=2-nsmet=5',cache=1)
    df = occ.df
    df = df[df.smetc.between(-0.4,0.4)]
    xk = 'smetc'
    dist = 'warm'

    size = 'se'
    cut = df.ix[dist,size]
    plot_rates(xk, cut, size)
    sampler = SamplerSmet('fit_smet-{}-{}'.format(dist,size),size)
    sampler.plot_all()
    
    size = 'sn'
    cut = df.ix[dist,size]
    plot_rates(xk, cut, size)
    sampler = SamplerSmet('fit_smet-{}-{}'.format(dist,size),size)
    sampler.plot_all()
    
    size = 'ss'
    cut = df.ix[dist,size]
    plot_rates(xk, cut, size)
    sampler = SamplerSmet('fit_smet-{}-{}'.format(dist,size),size)
    sampler.plot_all()
    
    size = 'jup'
    cut = df.ix[dist,size]
    plot_rates(xk, cut, size)
    sampler = SamplerSmet('fit_smet-{}-{}'.format(dist,size),size)
    #sampler.plot_all()

    ylabel('Planets Per 100 Stars')
    ylim(1e-3,1)


def fig_smet_hot():

    binwper = 0.25
    xk = 'smetc'
    key = 'occur-per={:f}-prad=physical-smet=0.2'.format(binwper)
    occ = cksmet.io.load_object(key)
    df = occ.df
    df = df[df.smetc.between(-0.4,0.4)]

    size = 'se'
    cut = df.query('1 < perc < 10 and 1.0 < pradc < 1.7')
    cut = sum_cells_per(cut)
    plot_rates(xk, cut, size)
    sampler = SamplerPerSmet('fit_persmet-hot-{}'.format(size),size)
    sampler.plot_all()

    cut = df.query('1 < perc < 10 and 1.7 < pradc < 4.0')
    cut = sum_cells_per(cut)
    size = 'sn'
    plot_rates(xk, cut, size)
    sampler = SamplerPerSmet('fit_persmet-hot-{}'.format(size),size)
    sampler.plot_all()

    occ = cksmet.io.load_object('occur-nper=2-nsmet=5',cache=1)
    df = occ.df
    df = df[df.smetc.between(-0.4,0.4)]

    size = 'ss'
    cut = df.ix['hot',size]
    plot_rates(xk, cut, size)
    sampler = SamplerPerSmet('fit_persmet-hot-{}'.format(size),size)
    sampler.plot_all()

    size = 'jup'
    cut = df.ix['hot',size]
    plot_rates(xk, cut, size)
    sampler = SamplerPerSmet('fit_persmet-hot-{}'.format(size),size)
    sampler.plot_all()
    

    ylabel('Planets Per 100 Stars')
    ylim(1e-5,1)


def yticks_planets_per_100_stars():
    yt = np.array([1e-4,3,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1,3e-1,1])
    syt = map(lambda x: x if x < 1 else "%.0f" % x , yt*100)
    yticks(yt,syt)

def fig_smet():
    fig, axL = subplots(ncols=2,nrows=1,figsize=(8,7),sharex=True,sharey=True)
    sca(axL[0])
    fig_smet_hot()

    yticks_planets_per_100_stars()
    title('$P = 1-10$ days')
    fig_label('a')

    sca(axL[1])
    fig_smet_warm()
    yticks_planets_per_100_stars()
    ylim(1e-4,1)
    ylabel('')
    fig_label('b')
    title('$P = 10-100$ days')
    
    fig.set_tight_layout(True)
    xlim(-0.4,0.4)
    setp(axL, xlabel='[Fe/H]')

    sca(axL[0])
    i = 1
    for size in 'se sn ss jup'.split():
        s =  "\n"*i + namedict[size] 
        text(0.1, 1, s, color=ptcolor[size],va='top',ha='left')
        i+=1

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
        bbox=dict(ec='none', fc='w', alpha=0.0)
    )
