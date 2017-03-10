from matplotlib.pylab import *
import seaborn as sns
import pandas as pd

import cksphys.io
import cksmet.io
import cksmet.cuts
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

import string
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
    fig,axL = subplots(nrows=3, ncols=2,figsize=(8,8))
    
    lamo = cksmet.io.load_table('lamost-dr2-cuts',cache=1)
    lamo = lamo.groupby('id_kic').first()
    lamoc = lamo.query('isany==False')

    cks = cksmet.io.load_table('cks-cuts',cache=1)
    cksc = cks.query('isany==False')
    
    cks = cks.groupby('id_kic').first()
    cksc = cksc.groupby('id_kic').first()

    plotkw = dict(marker='o',ms=1.5,lw=0)
    bins = np.arange(9,16.0001,0.2)
    histkw = dict(bins=bins,histtype='stepfilled')
    
    bins = np.arange(-1.0,0.5001,0.1)
    smethistkw = dict(bins=bins,histtype='stepfilled')

    sca(axL[0,0])
    plot(cks.cks_steff, cks.cks_slogg, color='LightGray',**plotkw)
    plot(cksc.cks_steff, cksc.cks_slogg,color='b',**plotkw)
    add_anchored('CKS',loc=1)

    sca(axL[0,1])
    plot(lamo.lamo_steff, lamo.lamo_slogg,color='LightGray',**plotkw)
    plot(lamoc.lamo_steff, lamoc.lamo_slogg,color='b',**plotkw)
    add_anchored('LAMOST',loc=1)

    sca(axL[1,0])
    hist(cks.koi_kepmag,color='LightGray',**histkw)
    hist(cksc.koi_kepmag,color='b',**histkw)
    add_anchored('CKS',loc=1)

    sca(axL[1,1])
    hist(lamo.koi_kepmag,color='LightGray',**histkw)
    hist(lamoc.koi_kepmag,color='b',**histkw)
    add_anchored('LAMOST',loc=1)
    
    sca(axL[2,0])
    hist(cks.cks_smet.dropna(),color='LightGray',**smethistkw)
    hist(cksc.cks_smet.dropna(),color='b',**smethistkw)
    add_anchored('CKS',loc=1)

    sca(axL[2,1])
    hist(lamo.lamo_smet,color='LightGray',**smethistkw)
    hist(lamoc.lamo_smet,color='b',**smethistkw)
    add_anchored('LAMOST',loc=1)

    setp(
        axL[0,:], xlim=(7000,4000), ylim=(5,3.0), xlabel='Teff (K)', 
        ylabel='logg (dex)'
    )
    fig.set_tight_layout(True)

    setp(axL[1,:],ylabel='Number of Stars',xlabel='Kepmag')
    setp(axL[2,:],ylabel='Number of Stars',xlabel='[Fe/H]')

    for ax, letter in zip(axL.flatten(),string.ascii_lowercase):
        sca(ax)
        add_anchored(
            letter,loc=2, frameon=False, 
            prop=dict(size='large', weight='bold')
        )
