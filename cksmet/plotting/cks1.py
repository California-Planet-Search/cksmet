from matplotlib.pylab import *
import seaborn as sns
import pandas as pd

import cksphys.io
import cksmet.io
import cksmet.cuts
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

figsize = (4,5)
errorbar_kw = dict(fmt='.',markersize=5,color='b')

sns.set_style('whitegrid')
sns.set_color_codes()

rpticks = [0.2, 0.3, 0.4, 0.5, 0.7, 1, 2, 3, 4, 5, 7, 10, 20]
texteff = '$\mathregular{T}_{\mathregular{eff}}$'
def add_anchored(*args,**kwargs):
    ax = gca()
    at = AnchoredText(*args,**kwargs)
    ax.add_artist(at)


def flicker():
    df = cksmet.io.load_table('cks+bast14')
    cksmet.plotting.cks1.comparison('slogg',df.cks_slogg,df.bast_slogg,'B14')
    pdffn = 'paper/fig_cks-bastien14-logg.pdf'
    print "saving {}".format(pdffn)
    gcf().savefig(pdffn)

    figure()
    pdffn = 'paper/fig_cks-bastien14-logg-kepmag.pdf'
    plot(df.koi_kepmag,df.bast_slogg-df.cks_slogg,'.')
    ylabel('$\Delta$ logg [CKS]')
    xlabel('Kepmag')

    print "saving {}".format(pdffn)
    gcf().savefig(pdffn)


def comparison(key,x1,x2,label1='CKS',label2='Comp',fig0=None, axL0=None):
    if key=='steff':
        x3 = x2 - x1 
        fig, axL = subplots_compare(x1, x2, x3, fig0=fig0,axL0=axL0,**errorbar_kw)
        sca(axL[0])
        ylabel('{} [{}] (K)'.format(texteff,label2))
        xlim(4500,7000)
        ylim(4500,7000)
        one2one()
        sca(axL[1])
        xlabel('{} [{}] (K)'.format(texteff,label1))
        ylabel('{} - {}'.format(label2,label1))
        ylim(-400,400)
        gca().yaxis.set_major_locator(MaxNLocator(5))
        axhline(0,color='g',zorder=1)

        sca(axL[0])
    if key=='slogg':
        x3 = x2 - x1 
        fig, axL = subplots_compare(x1, x2, x3, fig0=fig0,axL0=axL0,**errorbar_kw)
        sca(axL[0])
        ylabel('logg [{}] (dex)'.format(label2))
        xlim(3.5,5.0)
        ylim(3.5,5.0)
        one2one()

        sca(axL[1])
        xlabel('logg [{}] (dex)'.format(label1))
        ylabel('{} - {}'.format(label2,label1))
        ylim(-0.5,0.5)
        gca().yaxis.set_major_locator(MaxNLocator(5))
        axhline(0,color='g',zorder=1)

        sca(axL[0])

    if key=='smet':
        x3 = x2 - x1 
        fig, axL = subplots_compare(x1, x2, x3, fig0=fig0,axL0=axL0,**errorbar_kw)
        sca(axL[0])
        ylabel('[Fe/H] [{}] (dex)'.format(label2))
        xlim(-0.6,0.6)
        ylim(-0.6,0.6)
        one2one()

        sca(axL[1])
        xlabel('[Fe/H] [{}] (dex)'.format(label1))
        ylabel('{} - {}'.format(label2,label1))
        ylim(-0.3,0.3)
        gca().yaxis.set_major_locator(MaxNLocator(5))
        axhline(0,color='g',zorder=1)

        sca(axL[0])

    fmt = dict(param=key, mdiff=np.mean(x3),sdiff=np.std(x3),ndiff=len(x3))

    if key=='steff':
        s = "Mean($\Delta$) = {mdiff:.0f} K\nRMS($\Delta$) = {sdiff:.0f} K".format(
            **fmt
        )

    if key=='slogg':
        s = "Mean($\Delta$) = {mdiff:.2f} dex\nRMS($\Delta$) = {sdiff:.2f} dex".format(
            **fmt
        )

    if key=='smet':
        s = "Mean($\Delta$) = {mdiff:.3f} dex\nRMS($\Delta$) = {sdiff:.3f} dex".format(
            **fmt
        )
        
    print " ".join("{}: {}".format(k,fmt[k]) for k in 'param ndiff mdiff sdiff'.split())
    add_anchored(s,loc=2)

def comparison_three(table):
    fig = figure(figsize=(10,4.5))
    ax1 = subplot2grid((10,3), (0,0), rowspan=7)
    ax2 = subplot2grid((10,3), (7,0), rowspan=3, sharex=ax1)
    ax3 = subplot2grid((10,3), (0,1), rowspan=7)
    ax4 = subplot2grid((10,3), (7,1), rowspan=3, sharex=ax3)
    ax5 = subplot2grid((10,3), (0,2), rowspan=7)
    ax6 = subplot2grid((10,3), (7,2), rowspan=3, sharex=ax5)

    if table=='cks-buch14':
        df = cksmet.io.load_table('buch14-stars+cks')
        df = df.groupby('id_koi',as_index=False).last()
        steff1 = df.cks_steff
        steff2 = df.spc_steff
        slogg1 = df.cks_slogg
        slogg2 = df.spc_slogg
        smet1 = df.cks_smet
        smet2 = df.spc_smet
        kw = dict(label1='CKS',label2='B14',fig0=fig)
        comparison('steff',steff1,steff2, axL0=[ax1,ax2],**kw)
        comparison('slogg',slogg1,slogg2, axL0=[ax3,ax4],**kw)
        comparison('smet',smet1,smet2, axL0=[ax5,ax6],**kw)

    if table=='cks-lamost':
        df = cksmet.io.load_table('lamost-dr2+cks')
        df = df.groupby('id_koi',as_index=False).last()
        steff1 = df.cks_steff
        steff2 = df.lamo_steff
        slogg1 = df.cks_slogg
        slogg2 = df.lamo_slogg
        smet1 = df.cks_smet
        smet2 = df.lamo_smet
        kw = dict(label1='CKS',label2='LAMOST',fig0=fig)
        comparison('steff',steff1,steff2, axL0=[ax1,ax2],**kw)
        comparison('slogg',slogg1,slogg2, axL0=[ax3,ax4],**kw)
        comparison('smet',smet1,smet2, axL0=[ax5,ax6],**kw)

    if table=='cks-kic':
        df = cksmet.io.load_table('kic+cks')
        df = df.groupby('id_kic',as_index=False).last()
        df = df.query('kic_steff > 0')
        steff1 = df.cks_steff
        steff2 = df.kic_steff
        slogg1 = df.cks_slogg
        slogg2 = df.kic_slogg
        smet1 = df.cks_smet
        smet2 = df.kic_smet

        kw = dict(label1='CKS',label2='KIC',fig0=fig)
        comparison('steff',steff1,steff2, axL0=[ax1,ax2],**kw)
        comparison('slogg',slogg1,slogg2, axL0=[ax3,ax4],**kw)
        comparison('smet',smet1,smet2, axL0=[ax5,ax6],**kw)

    if table=='cks-huber14':
        df = cksmet.io.load_table('huber14+cks')
        df = df.groupby('id_kic',as_index=False).last()
        steff1 = df.cks_steff
        steff2 = df.huber_steff
        slogg1 = df.cks_slogg
        slogg2 = df.huber_slogg
        smet1 = df.cks_smet
        smet2 = df.huber_smet

        kw = dict(label1='CKS',label2='H14',fig0=fig)
        comparison('steff',steff1,steff2, axL0=[ax1,ax2],**kw)
        comparison('slogg',slogg1,slogg2, axL0=[ax3,ax4],**kw)
        comparison('smet',smet1,smet2, axL0=[ax5,ax6],**kw)

    if table=='cks-bruntt12':
        df = cksmet.io.load_table('bruntt12+cks')
        df = df.groupby('id_koi',as_index=False).last()
        steff1 = df.cks_steff
        steff2 = df.bruntt_steff
        slogg1 = df.cks_slogg
        slogg2 = df.bruntt_slogg
        smet1 = df.cks_smet
        smet2 = df.bruntt_smet

        kw = dict(label1='CKS',label2='B12',fig0=fig)
        comparison('steff',steff1,steff2, axL0=[ax1,ax2],**kw)
        comparison('slogg',slogg1,slogg2, axL0=[ax3,ax4],**kw)
        comparison('smet',smet1,smet2, axL0=[ax5,ax6],**kw)

    if table=='cks-everett13':
        df = cksmet.io.load_table('everett13+cks')
        df = df.groupby('id_koi',as_index=False).last()
        steff1 = df.cks_steff
        steff2 = df.everett_steff
        slogg1 = df.cks_slogg
        slogg2 = df.everett_slogg
        smet1 = df.cks_smet
        smet2 = df.everett_smet

        kw = dict(label1='CKS',label2='E13',fig0=fig)
        comparison('steff',steff1,steff2, axL0=[ax1,ax2],**kw)
        comparison('slogg',slogg1,slogg2, axL0=[ax3,ax4],**kw)
        comparison('smet',smet1,smet2, axL0=[ax5,ax6],**kw)

    if table=='cks-endl16':
        df = cksmet.io.load_table('endl16+cks')
        df = df.groupby('id_koi',as_index=False).last()
        steff1 = df.cks_steff
        steff2 = df.endl_steff
        slogg1 = df.cks_slogg
        slogg2 = df.endl_slogg
        smet1 = df.cks_smet
        smet2 = df.endl_smet

        kw = dict(label1='CKS',label2='E16',fig0=fig)
        comparison('steff',steff1,steff2, axL0=[ax1,ax2],**kw)
        comparison('slogg',slogg1,slogg2, axL0=[ax3,ax4],**kw)
        comparison('smet',smet1,smet2, axL0=[ax5,ax6],**kw)

    if table=='sm-bruntt12':
        df = cksmet.io.load_table('sm+bruntt12')
        df = df.groupby('id_kic',as_index=False).last()
        steff1 = df.sm_steff
        steff2 = df.bruntt_steff
        slogg1 = df.sm_slogg
        slogg2 = df.bruntt_slogg
        smet1 = df.sm_smet
        smet2 = df.bruntt_smet

        kw = dict(label1='SM',label2='B12',fig0=fig)
        comparison('steff',steff1,steff2, axL0=[ax1,ax2],**kw)
        comparison('slogg',slogg1,slogg2, axL0=[ax3,ax4],**kw)
        comparison('smet',smet1,smet2, axL0=[ax5,ax6],**kw)

    fig.set_tight_layout(True)
    
def cks_comparisons():
    plotnames = [
        'sm-bruntt12',
        'cks-endl16',
        'cks-bruntt12',
        'cks-endl16',
        'cks-everett13',
        'cks-huber14',
        'cks-kic',
        'cks-buch14',
        'cks-lamost',
    ]

    for plotname in plotnames:
        print "comparison: {}".format(plotname)
        comparison_three(plotname)
        gcf().savefig('paper/fig_{}.pdf'.format(plotname))

def provision_figure():
    fig = figure(figsize=figsize)
    ax1 = subplot2grid((4,1), (0,0), rowspan=3)
    ax2 = subplot2grid((4,1), (3,0), rowspan=1, sharex=ax1)
    axL = [ax1,ax2]
    return fig,axL 

def one2one(**kwargs):
    xl = xlim()
    plot(xl,xl,**kwargs)

def subplots_compare(x1, x2, x3, xerr3=None, fig0=None, axL0=None, **kwargs):
    if fig0 is None:
        fig, axL = provision_figure()
    else:
        fig = fig0
        axL = axL0

    fig.set_tight_layout(False)
    fig.subplots_adjust(hspace=0.4,left=0.17,top=0.95,right=0.90)
    sca(axL[0])
    grid()
    errorbar(x1,x2,**kwargs)
    sca(axL[1])
    grid() 
    errorbar(x1,x3,**kwargs)
    return fig,axL

