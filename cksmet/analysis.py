# Load up the culled CKS sample.
import cksmet.occur2
import cksmet.occur
import cksmet.io
import pandas as pd 

def compute_binned_occurrence(smet_bins=[-1.0, 0, 0.5]):
    print "Initializing completeness object"
    field = cksmet.io.load_table('field-cuts',cache=1)
    field = field.query('~isany')
    field = field.rename(columns={'huber_srad':'srad','huber_smass':'smass'})
    field.index = field.id_kic
    per_bins = cksmet.grid.per_bins_dict['xfine']
    prad_bins = cksmet.grid.prad_bins_dict['xfine']
    bins_dict = {'per': per_bins,'prad': prad_bins}
    spacing_dict = {'per':'log','prad':'log'}
    grid = cksmet.grid.Grid(bins_dict,spacing_dict)

    prob_det_mes_name = 'step-%i' % 12
    impact_transit = 0.9
    comp = cksmet.occur.Completeness(
        field, grid, prob_det_mes_name, impact_transit 
    )
    comp.mesfac = 1
    comp.compute_grid_prob_det()
    comp.compute_grid_prob_tr()
    comp.init_prob_det_interp()

    print "Initializing occurrence object"
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
    occur = cksmet.occur2.Occurrence(plnt, comp, nstars, smet_field=smet_field)

    print "Define occurrence grid"
    per_bins = cksmet.grid.per_bins_dict['four-per-decade']
    prad_bins = [1.0, 1.7, 4.0, 8.0, 24.0]

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
        s+="rate={rate:.4f} "
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
    return out
