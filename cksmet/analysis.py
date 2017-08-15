import numpy as np
import pandas as pd 

import cksmet.io
import cksmet.grid
import cksmet.comp
import cksmet.occur
import cksmet.fit
# Commonly used qualitative indecies
dist = pd.DataFrame(
    index=[1.0,10,300], data=['hot','warm','cool'], columns=['dist']
)
size = pd.DataFrame(
    index=[1.0,1.7,4.0,8.0],data=['se','sn','ss','jup',],columns=['size']
)
smet = pd.DataFrame(index=[-1,0],data=['sub','sup'],columns=['smet'])

# Some nice, convenient grids
per_bins_dict = {
# sqrt2 
    'xfine': np.round(np.logspace(np.log10(0.1),np.log10(1000),33),4),
# sqrt2 
    'fine': [ 
        1.00, 1.41, 2.00,  2.83,  4.00, 5.66, 8.00,  
        11.3, 16., 22.6, 32.0, 45.3, 64.0, 90.5, 128., 
        181,  256 ],
    'coarse': [ 
        1.00, 2.00,  4.00, 8.00,  
        16., 32.0, 64.0, 128., 256 ],
    'four-per-decade': 10**np.arange(-0.5,3.001,0.25)    
}    


prad_bins_dict = {
    'xfine': np.round(np.logspace(np.log10(0.5),np.log10(32),49 ),2),
    'fine': np.round(np.logspace(np.log10(0.5),np.log10(32),25 ),2),
    'two-per-octave': 2**np.arange(-1,5.001,0.5),
    'coarse': [ 
         0.5,  0.71, 1.00, 1.41, 2.0, 2.83, 4.00, 5.66, 8.0, 11.31, 16.0
    ]
}

def load_completeness():
    method = 'fulton-gamma' # treatment for planet detectability
    impact = 0.9 # maximum impact parameter considered.

    field = cksmet.io.load_table('field-cuts',cache=1)
    field = field.query('~isany')
    field = field.rename(columns={'huber_srad':'srad','huber_smass':'smass'})
    field.index = field.id_kic

    comp_per_bins = per_bins_dict['xfine']
    comp_prad_bins = prad_bins_dict['xfine']
    comp_bins_dict = {'per': comp_per_bins,'prad': comp_prad_bins}
    spacing_dict = {'per':'log','prad':'log'}
    grid = cksmet.grid.Grid(comp_bins_dict,spacing_dict)

    comp = cksmet.comp.Completeness(field, grid, method, impact)
    comp.compute_grid_prob_det()
    comp.compute_grid_prob_tr()
    comp.init_prob_det_interp()
    return comp

def compute_binned_occurrence(per_bins, prad_bins, smet_bins):
    print "Initializing occurrence object"
    comp = cksmet.io.load_object('comp',cache=1)
    cks = cksmet.io.load_table('cks-cuts')
    cks = cks[~cks.isany]
    cks = cks.dropna(subset=['iso_prad'])
    namemap = {
        'iso_srad':'srad','koi_period':'per','iso_prad':'prad',
        'iso_sma':'smax','cks_smet':'smet','koi_impact':'impact',
        'koi_max_mult_ev':'mes'
    }
    plnt = cks.rename(columns=namemap)

    lamo = cksmet.io.load_table('lamost-dr2-cal-cuts')
    lamo = lamo[~lamo.isany]
    smet_field = lamo.lamo_smet
    nstars = 35369.0
    occur = cksmet.occur.Occurrence(plnt, comp, nstars, smet_field=smet_field)

    print "Define occurrence grid"

    bins_dict = {'per': per_bins,'prad': prad_bins,'smet':smet_bins}
    spacing_dict = {'per':'log','prad':'log','smet':'linear'}
    grid = cksmet.grid.Grid(bins_dict,spacing_dict)

    print "Loop over occurrence grid"
    out = []
    i_count = 0
    def _print_row(d):
        s = ""
        s+="per={per1:.1f}-{per2:.1f}, "
        s+="prad={prad1:.1f}-{prad2:.1f}, "
        s+="smet={smet1:.1f}-{smet2:.1f}, "
        s+="rate={rate_str:s} "
        s = s.format(**_out)
        print s 

    for i, row in grid.ds.to_dataframe().iterrows():
        row = dict(row)
        _out = occur.occurence_box(row)
        _out = dict(row, **_out)
        out.append(_out)
        i_count+=1
        if i_count%10==0:
            _print_row(out)

    out = pd.DataFrame(out)
    for k in 'per1 perc per2 prad1 pradc prad2 smet1 smetc smet2'.split():
        out[k+'r'] = out[k].round(decimals=1)

    occur.df = out
    return occur

def load_occur(key):
    if key=='occur-test':
        per_bins = [1,10,100]
        prad_bins = [1,4,16]
        smet_bins = [-1, 0, 0.5]
        occ = compute_binned_occurrence(
            per_bins, prad_bins, smet_bins
        )

    elif key=='occur-nsmet=5':
        per_bins = per_bins_dict['four-per-decade']
        prad_bins = [1.0, 1.7, 4.0, 8.0, 24.0]
        smet_bins = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4] 
        occ = compute_binned_occurrence(per_bins, prad_bins, smet_bins)

        occ.df = pd.merge(occ.df,dist,left_on='per1',right_index=True)
        occ.df = pd.merge(occ.df,size,left_on='prad1',right_index=True)
        occ.df.set_index(['dist', 'size'], inplace=True)

    elif key=='occur-nsmet=2':
        per_bins = per_bins_dict['four-per-decade']
        prad_bins = [1.0, 1.7, 4.0, 8.0, 24.0]
        smet_bins = [-1, 0, 0.5] 
        occ = compute_binned_occurrence(
            per_bins, prad_bins, smet_bins
        )
        occ.df = pd.merge(occ.df, size, left_on='prad1',right_index=True)
        occ.df = pd.merge(occ.df, smet, left_on='smet1',right_index=True)
        occ.df.set_index(['size', 'smet'], inplace=True)

    elif key=='occur-nper=2-nsmet=5':
        per_bins = [1, 10, 100]
        prad_bins = [1.0, 1.7, 4.0, 8.0, 24.0]
        smet_bins = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4] 
        occ = compute_binned_occurrence(
            per_bins, prad_bins, smet_bins
        )
        occ.df = pd.merge(occ.df,dist,left_on='per1',right_index=True)
        occ.df = pd.merge(occ.df,size,left_on='prad1',right_index=True)
        occ.df.set_index(['dist', 'size'], inplace=True)

    elif key=='occur-nsmet=1':
        per_bins = per_bins_dict['four-per-decade']
        prad_bins = prad_bins_dict['two-per-octave']
        smet_bins = [-1, 0.5]
        occ = compute_binned_occurrence(
            per_bins, prad_bins, smet_bins
        )
    
    elif key=='occur-hot-jup':
        per_bins = [1,10]
        prad_bins = [1,24]
        smet_bins = [-1, 0.5]
        occ = compute_binned_occurrence(
            per_bins, prad_bins, smet_bins
        )
    
    return occ

def load_fit(key):
    if key.count('fit_smet')==1:
        # Fit hot SN
        _, dist, size = key.split('-')
        df = cksmet.io.load_table('occur-nper=2-nsmet=5',cache=1)
        cut = df.ix[dist,size]
        cut = cut[cut.smetc.between(-0.4,0.4)]
        fit = cksmet.fit.Exponential(cut.smetc, cut.nplnt, cut.ntrial)
        fit.fit()
        fit.mcmc(burn=300, steps=600, thin=1, nwalkers=100)
        fit.print_parameters()

    elif key.count('fit_per')==1:
        # Fit hot SN

        _, smet, size = key.split('-')
        df = cksmet.io.load_table('occur-nsmet=2',cache=1)
        cut = df.ix[size,smet]
        fit = cksmet.fit.PowerLawCutoff(cut.perc, cut.nplnt, cut.ntrial)
        cut = cut[cut.perc.between(1,350) & (cut.prob_det_mean > 0.25)]

        fit.fit()
        fit.mcmc(burn=300, steps=600, thin=1, nwalkers=100)
        fit.print_parameters()

    else:
        assert False, "{} not supported".format(key)

    fit.occur = cut
    return fit
