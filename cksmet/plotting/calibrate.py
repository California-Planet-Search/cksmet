import seaborn as sns

import cksmet.io
import cksspec.plotting.catalogs
from cksspec.plotting.compare import comparison

import pandas as pd
from matplotlib.pylab import *

def load_comparison(table):
    if table=='lamost-dr2':
        df = cksmet.io.load_table(table,cache=1)
    else:
        df = cksmet.io.load_table(table)

    namemap = {
        'lamo_steff':'teff','lamo_slogg':'logg','lamo_smet':'fe'
    }
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


def validation_lamo():
    sns.set(
        style='ticks',
        rc={'ytick.major.size':3.0,
            'xtick.major.size':3.0,
            'xtick.direction': u'in',
            'ytick.direction': u'in',}
    )

    df = load_comparison('lamost-dr2')
    dfcal = load_comparison('lamost-dr2-cal')

    # Raw
    cksspec.plotting.catalogs.fivepane(
        df,suffixes=['_lib',''], color='fe',
        kw_twopane=dict(labels=['CKS','LAMO'])
    )
    plot_prefix = 'lamost-dr2-on-cks'
    gcf().get_axes()[0].set_title(plot_prefix+' (raw)')
    pdffn = 'fig_{}-raw.pdf'.format(plot_prefix)
    gcf().savefig(pdffn)
    print pdffn

    # Cal
    cksspec.plotting.catalogs.fivepane(
        dfcal,suffixes=['_lib',''], color='fe',
        kw_twopane=dict(labels=['CKS','LAMO'])
    )
    gcf().get_axes()[0].set_title(plot_prefix+' (cal)')
    pdffn = 'fig_{}-cal.pdf'.format(plot_prefix)
    gcf().savefig(pdffn)
    print pdffn

    fig = figure(figsize=(7,3.5))
    shape = (10,2)
    ax1 = subplot2grid(shape, (0,0), rowspan=8)
    ax2 = subplot2grid(shape, (8,0), rowspan=2, sharex=ax1)
    ax3 = subplot2grid(shape, (0,1), rowspan=8)
    ax4 = subplot2grid(shape, (8,1), rowspan=2, sharex=ax3)

    kw = dict(label1='CKS',label2='LAMOST',fig0=fig)
    comparison('smet',df.fe_lib,df.fe, axL0=[ax1,ax2],**kw)
    comparison('smet',dfcal.fe_lib,dfcal.fe, axL0=[ax3,ax4],**kw)

    fig.set_tight_layout(True)
    fig.set_tight_layout(False)
    fig.subplots_adjust(hspace=0.001,left=0.12,top=0.96,right=0.98,wspace=0.4,bottom=0.14)


