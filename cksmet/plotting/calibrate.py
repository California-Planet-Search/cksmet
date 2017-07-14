import seaborn as sns

import cksmet.io
import cksspec.plotting.catalogs
from cksspec.plotting.compare import comparison

import pandas as pd
from matplotlib.pylab import *

def validation_lamo():
    sns.set(style='ticks',rc={'ytick.major.size':3.0,'xtick.major.size':3.0,'xtick.direction': u'in','ytick.direction': u'in',})

    df = cksmet.io.load_table('lamost-dr2-cal')

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

    cksspec.plotting.catalogs.fivepane(
        df,suffixes=['_lib',''], color='fe',
        kw_twopane=dict(labels=['CKS','LAMO'])
    )
    plot_prefix = 'lamost-dr2-cal-on-cks-cal'
    gcf().get_axes()[0].set_title(plot_prefix+' (cal)')
    pdffn = 'fig_{}-cal.pdf'.format(plot_prefix)
    gcf().savefig(pdffn)
    print pdffn


    comparison('smet',df.fe_lib,df.fe,label1='CKS',label2='LAMO')
    fig = gcf()
    fig.set_tight_layout(True)
    draw()
    fig.set_tight_layout(False)
    fig.subplots_adjust(hspace=0.001,left=0.2,top=0.96,right=0.96,bottom=0.15)
    fig.savefig('fig_sm-lamo.pdf')

