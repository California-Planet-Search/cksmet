import string

from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib.pylab import *
import seaborn as sns
import pandas as pd

import cksmet.io
import cksmet.cuts
from astropy import constants as c

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

def samples():
    ncols = 3
    nrows = 3
    height = 2.5 * nrows
    width = 3 * ncols
    fig, axL = subplots(nrows=nrows, ncols=ncols, figsize=(width, height))

    lamo = cksmet.io.load_table('lamost-dr2-cal-cuts',cache=2)
    lamoc = lamo[~lamo.isany]

    cks = cksmet.io.load_table('cks-cuts',cache=1)
    cksc = cks.query('isany==False')
    
    cks = cks.groupby('id_kic').first()
    cksc = cksc.groupby('id_kic').first()

    
    field = cksmet.io.load_table('huber14+cdpp',cache=1)

    field = cksmet.cuts.add_cuts(field, cksmet.cuts.plnt_cuttypes, 'field')
    fieldc = field.query('isany==False')

    plotkw = dict(marker='o',ms=1.5,lw=0)
    bins = np.arange(9,16.0001,0.2)
    histkw = dict(bins=bins,histtype='stepfilled')
    
    bins = np.arange(-1.0,0.5001,0.1)
    smethistkw = dict(bins=bins,histtype='stepfilled')

    jcks = 0
    jfield = 1 
    jlamo = 2

    # HR diagrams
    sca(axL[0,jcks])
    plot(cks.cks_steff, cks.cks_slogg, color='LightGray',**plotkw)
    plot(cksc.cks_steff, cksc.cks_slogg,color='b',**plotkw)

    sca(axL[0,jfield])
    print list(field.columns)
    plot(field.huber_steff, field.huber_slogg, color='LightGray',rasterized=True,**plotkw)
    plot(fieldc.huber_steff, fieldc.huber_slogg,color='b',rasterized=True,**plotkw)

    sca(axL[0,2])
    plot(lamo.lamo_steff, lamo.lamo_slogg,color='LightGray',**plotkw)
    plot(lamoc.lamo_steff, lamoc.lamo_slogg,color='b',**plotkw)

    # KepMag
    sca(axL[1,jcks])
    hist(cks.kic_kepmag,color='LightGray',**histkw)
    hist(cksc.kic_kepmag,color='b',**histkw)

    sca(axL[1,jfield])
    hist(field.kepmag,color='LightGray',**histkw)
    hist(fieldc.kepmag,color='b',**histkw)

    sca(axL[1,2])
    hist(lamo.kic_kepmag,color='LightGray',**histkw)
    hist(lamoc.kic_kepmag,color='b',**histkw)
    
    # Metal
    
    sca(axL[2,jcks])
    hist(cks.cks_smet.dropna(),color='LightGray',**smethistkw)
    hist(cksc.cks_smet.dropna(),color='b',**smethistkw)


    sca(axL[2,2])
    hist(lamo.lamo_smet,color='LightGray',**smethistkw)
    hist(lamoc.lamo_smet,color='b',**smethistkw)

    setp(
        axL[0,:], xlim=(7000,4000), ylim=(5,3.0), xlabel='Teff (K)', 
        ylabel='logg (dex)'
    )
    fig.set_tight_layout(True)

    setp(axL[1,:],ylabel='Number of Stars',xlabel='Kepmag')
    setp(axL[2,:],ylabel='Number of Stars',xlabel='[Fe/H]')

    for ax in axL[:,2]:
        sca(ax)
        add_anchored('LAMOST',loc=1)
        
    for ax in axL[:,jcks]:
        sca(ax)
        add_anchored('CKS',loc=1)

    for ax in axL[:,jfield]:
        sca(ax)
        add_anchored('Field',loc=1)


    for ax, letter in zip(axL.flatten(),string.ascii_lowercase):
        sca(ax)
        add_anchored(
            letter,loc=2, frameon=False, 
            prop=dict(size='large', weight='bold')
        )

def lamo_detectability():
    querys = ['lamo_smet < 0','0 < lamo_smet']
    colors = ['b','r']
    labels = ['[Fe/H] < 0','[Fe/H] > 0']

    lamo = cksmet.io.load_table('lamost-dr2-cal-cuts+cdpp',cache=1)
    lamo = lamo[~lamo.isany]
    fig1 = figure(figsize=(4,4))
    semilogy()

    bins = arange(11,14.1,0.5)
    #querys = ['-0.5 < lamo_smet < -0.1','0.1 < lamo_smet < 0.5']
    for i in range(len(querys)):
        query = querys[i]
        color=colors[i]
        label=labels[i]
        cut = lamo.query(query) 
        plot(cut.kepmag,cut.cdpp3,'.',color=color,zorder=1,label=label,alpha=0.8)
        g = cut.groupby(pd.cut(cut.kepmag,bins))
        i = 0

        for interval, cdpp in g['cdpp3'].median().iteritems():
            if i==0:
                label=query
            else:
                label=None
            plot([interval.left,interval.right], [cdpp,cdpp],color='w',lw=6)
            plot([interval.left,interval.right], [cdpp,cdpp],color=color,lw=3)
            i+=1
            print query, interval.left,interval.right, cdpp

    xlim(11,14)
    ylim(15,150)
    legend()

    xlabel('kepmag')
    ylabel('CDPP3 (ppm)')
    minorticks_off()
    yticks([20,30,40,50,60,70,80,90,100],[20,30,40,50,60,70,80,90,100])
    fig1.set_tight_layout(True)


    fig2 = figure(figsize=(4,4))
    kepmag = (12,14)
    cut2 = lamo[lamo.kepmag.between(*kepmag)]
    fsamp =  100 * len(cut2) / len(lamo)
    print "for kepmag = {}, {:.2f}\% of stellar sample median CDPP3 is ".format(kepmag, fsamp)

    fig2 = figure(figsize=(4,4))
    for i in range(2):
        query = querys[i]
        color=colors[i]
        label=labels[i]
        cut = lamo.query(query) 
        plot(cut.lamo_steff,cut.lamo_slogg,'.',color=color)

    xlim(6500,4700)
    ylim(5.0,4.0)
    xlabel('Effective Temp. (K)')
    ylabel('logg (dex)')
    fig2.set_tight_layout(True)
    return fig1,fig2


def smet_snr():
    df = pd.read_csv('isoclassify-lamost-dr2.csv')
    huber14 = cksmet.io.load_table('huber14+cdpp',cache=1)
    lamo = cksmet.io.load_table('lamost-dr2-cal-cuts',cache=1)
    lamo = lamo[~lamo.isany]
    lamo = pd.merge(lamo,huber14['id_kic kepmag cdpp3'.split()],on='id_kic')
    df = pd.merge(lamo,df,left_on='id_kic',right_on='id_starname')

    fig = figure(figsize=(4,4))
    semilogy()
    srad = df.iso_srad * c.R_sun
    prad = 1 * c.R_earth
    df['snr'] = (prad / srad)**2 / (1e-6 * df.cdpp3)
    plot(df.lamo_smet,df.snr,',')

    bins = arange(-0.3,0.31,0.15,)
    g = df.groupby(pd.cut(df.lamo_smet,bins))
    i = 0
    for interval, snr in g['snr'].median().iteritems():
        x = [interval.left,interval.right]
        y = [snr,snr]
        plot(x,y,color='w',lw=6)
        plot(x,y,color='b',lw=3)
        i+=1
        print interval.left,interval.right, snr

    xlim(-0.4,0.4)
    ylim(0.1,10)
    yt = [0.1,0.3,1,3,10]
    yticks(yt,yt)
    xlabel('[Fe/H] (dex)')
    ylabel('Single Transit SNR')
    fig.set_tight_layout(True)
    return fig
