import pandas as pd 

import cksmet.fit
import cksmet.occur
import cksmet.io
import cksmet.grid
import cksmet.comp

def calc_completeness():
    field = cksmet.io.load_table('field-cuts',cache=1)
    field = field.query('~isany')
    field = field.rename(columns={'huber_srad':'srad','huber_smass':'smass'})
    field.index = field.id_kic

    comp_per_bins = cksmet.grid.per_bins_dict['xfine']
    comp_prad_bins = cksmet.grid.prad_bins_dict['xfine']
    comp_bins_dict = {'per': comp_per_bins,'prad': comp_prad_bins}
    spacing_dict = {'per':'log','prad':'log'}
    grid = cksmet.grid.Grid(comp_bins_dict,spacing_dict)

    prob_det_mes_name = 'step-%i' % 12
    impact_transit = 0.9
    comp = cksmet.comp.Completeness(
        field, grid, prob_det_mes_name, impact_transit 
    )
    comp.mesfac = 1
    comp.compute_grid_prob_det()
    comp.compute_grid_prob_tr()
    comp.init_prob_det_interp()
    return comp

def compute_binned_occurrence(per_bins, prad_bins, smet_bins):
    print "Initializing occurrence object"
    comp = cksmet.io.load_comp()
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

    return out

def fit_occurrence(key):
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
