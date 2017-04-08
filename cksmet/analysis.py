import cksmet.occur
import cksmet.io
from scipy.interpolate import interp1d
import cksmet.bins
import numpy as np
lamo = cksmet.io.load_table('lamost-dr2-cuts',cache=1)
lamo = lamo[~lamo.isany]
x = lamo.lamo_feh
xs = np.sort(x)
ys = np.arange(1, len(xs)+1)/float(len(xs))
metcdf = interp1d(xs,ys)
namemap = {
    'iso_srad':'srad','koi_period':'per','iso_prad':'prad','iso_sma':'smax',
    'koi_impact':'impact','koi_max_mult_ev':'mes'
}

def load_occur(smetbin=None):
    """
    Load occurrence object

    Args:
        smet (Optional[list]): lower and upper bounds on metallicity
    
    """

    cks = cksmet.io.load_table('cks-cuts')
    cks = cks[~cks.isany | cks.isgrazing]
    cks = cks.dropna(subset=['iso_prad'])


    # Constants
    mes = 12
    impact_transit = 0.8

    _metcdf = metcdf(smetbin)
    fstars = _metcdf[1] - _metcdf[0]
    nstars0 = 34654 
    nstars = nstars0 * fstars
    print "{}/{} stars ({:.3f} of sample) has fe = {} - {}".format(
        nstars, nstars0, fstars, *smetbin
    )


    idxfe = cks[cks.cks_smet.between(*smetbin)].index
    print "{}/{} hosts with smet = {} - {}".format(len(idxfe),len(cks),*smetbin)
    cks = cks.ix[idxfe]

    '''
    per_bins = cksmet.occur.per_bins_dict['xfine']
    prad_bins = cksmet.occur.prad_bins_dict['xfine']
    '''

    per_bins = cksmet.bins.per_bins_dict['xfine-hj']
    prad_bins = cksmet.bins.prad_bins_dict['xfine-hj']
    
    bins_dict = {'per': per_bins,'prad': prad_bins}
    print bins_dict 
    spacing_dict = {'per':'log','prad':'log'}
    grid = cksmet.occur.Grid(bins_dict,spacing_dict)
    
    field = cksmet.io.load_table('field-cuts',cache=1)
    field = field.query('~isany')
    field = field.rename(columns={'huber_srad':'srad','huber_smass':'smass'})
    #field['srad'] *= 1.0
    field.index = field.id_kic
    prob_det_mes_name = 'step-%i' % mes
    print "prob_det_mes_name {}".format(prob_det_mes_name)

    comp = cksmet.occur.Completeness(
        field, grid, prob_det_mes_name, impact_transit
    )
    comp.mesfac = 0.7469
    plnt = cks.rename(columns=namemap)
    
    #comp.compute_mes_factor(plnt.iloc[::10])
    comp.compute_grid_prob_det()
    comp.init_prob_det_interp()
    occ = cksmet.occur.Occurrence(plnt, nstars ,comp, grid)
    return occ



# Scan occurrence

'''
namemap = {
    'iso_srad':'srad','koi_period':'per','iso_prad':'prad','iso_sma':'smax',
    'koi_impact':'impact','koi_max_mult_ev':'mes'
}

field = cksmet.io.load_table('field-cuts',cache=1)
field = field.query('~isany')
field = field.rename(columns={'huber_srad':'srad','huber_smass':'smass'})
field['srad'] *= 1.0

field.index = field.id_kic

nstars = 34654
impact_transit = 0.7

for mes in [10,12,15,20]:
    prob_det_mes_name = 'step-%i' % mes
    print "{}".format(prob_det_mes_name)

    comp = cksmet.occur.Completeness(field, grid, prob_det_mes_name ,impact_transit)
    plnt = cks.rename(columns=namemap)
    #comp.compute_mes_factor(plnt.iloc[::10])
    comp.mesfac = 0.7469
    comp.compute_grid_prob_det()
    comp.init_prob_det_interp()
    occ = cksmet.occur.Occurrence(plnt, nstars ,comp, grid)
    ds = occ.compute_occurrence()

    # Checking
    dscut = ds.where((2 <= ds.prad1) & (ds.prad2 <= 4) & (10 <= ds.per1 ) & (ds.per2 <= 30))
    fcell = float(dscut.plnt_occur.sum(dim=['per','prad']))
    print "P = 1-30d, Rp = 2-4 Re  P17: {:.3f} (PPS)".format(fcell)

    # Hot juptiers
    dscut = ds.where((8 <= ds.prad1) & (ds.prad2 <= 32) & (1 <= ds.per1 ) & (ds.per2 <= 10))
    fcell = float(dscut.plnt_occur.sum(dim=['per','prad']))
    print "P = 1-10d, Rp = 8-32 Re  P17: {:.3f}, H12 0.4 % (PPS)".format(fcell)
    
    # H12 
    dscut = ds.where((2 <= ds.prad1) & (ds.prad2 <= 4) & (1 <= ds.per1 ) & (ds.per2 <= 50))
    fcell = float(dscut.plnt_occur.sum(dim=['per','prad']))
    print "P = 1-50d, Rp = 2-4 Re  P17: {:.3f} (PPS), H12: {:.3f} (PPS)".format(fcell,0.14)

    # F13
    dscut = ds.where((2 <= ds.prad1) & (ds.prad2 <= 4) & (1 <= ds.per1 ) & (ds.per2 <= 85))
    fcell = float(dscut.plnt_occur.sum(dim=['per','prad']))
    print "P = 1-85d, Rp = 2-4 Re, P17: {:.3f} (PPS), F13: {:.3f} (PPS) ".format(fcell,.23)

    # F13
    dscut = ds.where((1 <= ds.prad1) & (ds.prad2 <= 2) & (1 <= ds.per1 ) & (ds.per2 <= 85))
    fcell = float(dscut.plnt_occur.sum(dim=['per','prad']))
    print "P = 1-85d, Rp = 1-2 Re, P17: {:.3f} (PPS), F13: {:.3f} (PPS)".format(fcell,.33)

    # P13
    dscut = ds.where((2 <= ds.prad1) & (ds.prad2 <= 4) & (5 <= ds.per1 ) & (ds.per2 <= 100))
    fcell = float(dscut.plnt_occur.sum(dim=['per','prad']))
    print "P = 1-100d, Rp = 2-4 Re, P17: {:.3f} (PPS), P13: {:.3f} (FSP)".format(fcell,.245)

    # P13
    dscut = ds.where((1 <= ds.prad1) & (ds.prad2 <= 2) & (5 <= ds.per1 ) & (ds.per2 <= 100))
    fcell = float(dscut.plnt_occur.sum(dim=['per','prad']))
    print "P = 1-100d, Rp = 1-2 Re, P17: {:.3f} (PPS), P13: {:.3f} (FSP)".format(fcell,.262)

    # H12
    dscut = ds.where((1 <= ds.prad1) & (ds.prad2 <= 2) & (5 <= ds.per1 ) & (ds.per2 <= 100))
    fcell = float(dscut.plnt_occur.sum(dim=['per','prad']))
    print "P = 1-100d, Rp = 1-2 Re, P17: {:.3f} (PPS), P13: {:.3f} (FSP)".format(fcell,.262)
'''
