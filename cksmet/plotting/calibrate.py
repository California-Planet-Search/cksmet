import seaborn as sns

import cksmet.io
import cksspec.plotting.catalogs
from cksspec.plotting.compare import comparison

import pandas as pd
from matplotlib.pylab import *
from cksmet.plotting.config import *

def load_comparison(table):
    if table=='lamost-dr2':
        df = cksmet.io.load_table(table,cache=1)
    else:
        df = cksmet.io.load_table(table)

    namemap = {'lamo_steff':'teff','lamo_slogg':'logg','lamo_smet':'fe'}
    df = df.rename(columns=namemap)

    # Make fivepane plots
    lib = cksmet.io.load_table('cks')
    namemap = {
        'cks_steff':'teff_lib','cks_slogg':'logg_lib','cks_smet':'fe_lib'
    }
    lib = lib.rename(columns=namemap)
    df = pd.merge(df, lib, on=['id_kic'])
    df['teff_diff'] = df['teff'] - df['teff_lib'] 
    df['logg_diff'] = df['logg'] - df['logg_lib'] 
    df['fe_diff'] = df['fe'] - df['fe_lib'] 
    return df

def lamo_on_cks():
    sns.set(
        style='ticks',
        rc={'ytick.major.size':3.0,
            'xtick.major.size':3.0,
            'xtick.direction': u'in',
            'ytick.direction': u'in',}
    )
    sns.set_context('paper',font_scale=1.1)
    df = cksmet.io.load_table('lamost-cks-calibration-sample')

    dfcal = df.copy()
    calfn  = 'cal_lamo-to-cks.fits'
    namemap = {'teff_new':'teff','logg_new':'logg','fe_new':'fe'}
    dfcal = dfcal.rename(columns=namemap)

    dfcal = cksmet._calibrate.calibrate(dfcal, calfn, mode='uncal')
    dfcal = dfcal.drop(['delta'],axis=1)

    fig = figure(figsize=(7,3.5))
    shape = (10,2)
    ax1 = subplot2grid(shape, (0,0), rowspan=8)
    ax2 = subplot2grid(shape, (8,0), rowspan=2, sharex=ax1)
    ax3 = subplot2grid(shape, (0,1), rowspan=8)
    ax4 = subplot2grid(shape, (8,1), rowspan=2, sharex=ax3)

    kw = dict(label1='CKS',label2='LAMOST',fig0=fig)
    comparison('smet',df.fe_lib,df.fe_new, axL0=[ax1,ax2],**kw)
    comparison('smet',dfcal.fe_lib,dfcal.fe, axL0=[ax3,ax4],**kw)
    fig.set_tight_layout(True)
    fig.set_tight_layout(False)
    fig.subplots_adjust(
        hspace=0.001,left=0.12,top=0.96,right=0.98,wspace=0.4,bottom=0.14
    )

def fig_lamo_diff():
    import seaborn as sns
    sns.set(
        style='ticks',
        rc={'ytick.major.size':3.0,
            'xtick.major.size':3.0,
            'xtick.direction': u'in',
            'ytick.direction': u'in',}
    )
    sns.set_context('paper',font_scale=1.1)
    sns.set_style('ticks')
    df = cksmet.io.load_table('lamost-cks-calibration-sample-noclip')
    fig,axL = subplots(ncols=3,figsize=(8,3),sharey=True)
    axL= axL.flatten()

    kw = dict(lw=0,marker='.')
    sca(axL[0])
    fig_label("a")
    plot(df.fe_lib,df.fe_diff,**kw)
    xlabel('[Fe/H] (CKS)')
    grid()

    sca(axL[1])
    fig_label("b")
    plot(df.kic_kepmag,df.fe_diff,**kw)
    xlabel('$Kp$')
    grid()

    sca(axL[2])
    fig_label("c")
    plot(df.snrg,df.fe_diff,**kw)
    semilogx()
    xlim(400,8)
    xlabel('LAMOST (SNR)')
    xticks([300,100,30,10],[300,100,30,10])

    setp(axL,ylim=(-0.5,0.5))
    setp(axL[0],ylabel='$\Delta$ [Fe/H] (LAMOST-CKS)')
    tight_layout()
    grid()


    
