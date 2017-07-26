import cksmet.occur
import cksmet.io
from scipy.interpolate import interp1d
import numpy as np
import cksmet.grid
import cPickle as pickle
lamo = cksmet.io.load_table('lamost-dr2-cal-cuts',cache=1)
lamo = lamo[~lamo.isany]
x = lamo.lamo_smet
xs = np.sort(x)
ys = np.arange(1, len(xs)+1)/float(len(xs))
metcdf = interp1d(xs,ys)
namemap = {
    'iso_srad':'srad','koi_period':'per','iso_prad':'prad','iso_sma':'smax',
    'koi_impact':'impact','koi_max_mult_ev':'mes'
}

def cachefn(smetbin):
    return "occ_smet={:+.2f},{:+.2f}.pkl".format(*smetbin)


def load_occur(cachefn=None,cache=0,smetbin=None):
    """
    Load occurrence object

    Args:
        smet (Optional[list]): lower and upper bounds on metallicity
    
    """

    if cache==1:
        assert type(cachefn)!=type(None),"must set cachefn"

    if cache==1:
        with open(cachefn,'r') as f:
            occ = pickle.load(f)
            return occ

    # Load up the culled CKS sample.
    cks = cksmet.io.load_table('cks-cuts')
    cks = cks[~cks.isany]
    cks = cks.dropna(subset=['iso_prad'])

    # Constants
    mes = 12
    impact_transit = 0.9

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

    per_bins = cksmet.grid.per_bins_dict['xfine']
    prad_bins = cksmet.grid.prad_bins_dict['xfine']
    bins_dict = {'per': per_bins,'prad': prad_bins}
    spacing_dict = {'per':'log','prad':'log'}
    grid = cksmet.grid.Grid(bins_dict,spacing_dict)
    
    field = cksmet.io.load_table('field-cuts',cache=1)
    field = field.query('~isany')
    field = field.rename(columns={'huber_srad':'srad','huber_smass':'smass'})
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
    occ.set_grid_fine(grid)

    downsamp = {'per':2,'prad':4}
    occ.set_grid_loguni(downsamp)
    if cache==2:
        with open(cachefn,'w') as f:
            pickle.dump(occ,f)
            print "saving {} ".format(cachefn)
        
    return occ
