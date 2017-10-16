from matplotlib.pylab import *
import seaborn as sns
import matplotlib.patheffects as path_effects
from matplotlib.transforms import blended_transform_factory as btf
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import pandas as pd 

import cksmet.analysis
import cksmet.occur

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


def fig_checkerboard():
    occ  = cksmet.io.load_object('occur-nsmet=1') 

    prob_det_min = 0.1
    epsilon = 0.0001

    fig,ax = subplots(figsize=(7,8))

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
    cbar = colorbar(orientation='horizontal',ticks=[-4,-3,-2,-1])
    cbar.set_label('$\log_{10} (f_{cell})$')

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
    xlim(0.3,300)
    ylim(0.5, 32)
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
    _bdcolor = bdcolor[fmtkey]
    fill_between(peri, p16, p84, color=_bdcolor,zorder=1)
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
        bbox=dict(ec='none', fc='w', alpha=0.8)
    )
