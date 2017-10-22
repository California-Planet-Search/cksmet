import cksmet.stats
from statsmodels.distributions.empirical_distribution import ECDF

class Occurrence(object):
    def __init__(self, plnt, comp, nstars, smet_field=None):
        self.plnt = plnt
        self.comp = comp
        self.nstars = nstars
        self.plnt = plnt
        if smet_field is not None:
            self.smet_field = smet_field

    def occurence_box(self, limits):
        """Compute occurrence in a little box

        We make the assumption that the dN/dlogP and dN/lopRp is constant
        within a box.

        Args:
            limits (dict): must contain, per1, per2, prad1, prad2, can
                optionally contain, smet1, smet2
        """
        out = dict()
        prad1 = limits['prad1']
        prad2 = limits['prad2']
        per1 = limits['per1']
        per2 = limits['per2']

        # Get planet sample and number of stars
        cut = self.plnt.copy()
        cut = cut[cut.prad.between(prad1,prad2)]
        cut = cut[cut.per.between(per1,per2)]
        nstars = self.nstars
        if limits.has_key('smet1'):
            assert self.smet_field is not None, "Must set provide field met" 
            smet1 = limits['smet1']
            smet2 = limits['smet2']
            cut = cut[cut.smet.between(smet1,smet2)]
            ecdf = ECDF(self.smet_field)
            ecdf_smet1, ecdf_smet2 = ecdf([smet1,smet2])
            nstars = (ecdf_smet2 - ecdf_smet1) * nstars

        nplnt = len(cut)

        prob_trdet_mean, prob_det_mean = self.comp.mean_prob_trdet(
            per1, per2, prad1, prad2
        )
        ntrial = nstars * prob_trdet_mean
        rate = nplnt / ntrial

        nsample = int(1e4)
        binom = cksmet.stats.Binomial(ntrial, nplnt)
        samples = binom.sample(nsample) 

        uplim = nplnt==0
        rate = cksmet.stats.samples_to_rate(samples,uplim=uplim)
        out['ntrial'] = ntrial
        out['nplnt'] = nplnt
        out['prob_trdet_mean'] = prob_trdet_mean
        out['prob_det_mean'] = prob_det_mean
        out = dict(out,**rate)
        return out

