import numpy as np
from lmfit import Parameters 
import lmfit
import cksmet.model
import copy
fmt = {'kp':"{:.2f}",'beta':"{:.2f}",'gamma':"{:.1f}",'per0':"{:.1f}"}

class Fitter(object):
    def __init__(self, x, nplnt, ntrial):
        self.x = np.array(x)
        self.nplnt = np.array(nplnt)
        self.ntrial = np.array(ntrial)

    def _loglike(self, params):
        _loglike = cksmet.model.loglike_1d(
            params, self.model, self.x, self.nplnt, self.ntrial
        )
        return _loglike 

    def _negloglike(self,params):
        return -1.0 * self._loglike(params)

    def print_parameters(self):
        chain = self.flatchain
        q = chain.quantile([0.16,0.50,0.84])
        for k in q.columns:
            s = fmt[k]
            val = s.format(q.ix[0.50,k])
            err1 = s.format(q.ix[0.84,k] - q.ix[0.50,k])
            err2 = s.format(q.ix[0.16,k] - q.ix[0.50,k])
            s = "$%s^{+%s}_{%s}$" % (val, err1,err2)
            print r"{{{}}}{{{}}}".format(k,s)

    def sample(self, nsamp, x):
        psamples = self.flatchain.sample(nsamp)
        model_samples = []
        for i, row in psamples.iterrows():
            psamp = copy.copy(self.p0)
            for k in self.p0.valuesdict().keys():
                psamp[k].value = row[k]

            model_samples.append( self.model(psamp, x) ) 

        model_samples = np.vstack(model_samples)
        return model_samples

    def fit(self):
        res = lmfit.minimize(self._negloglike, self.p0, method='Nelder')
        lmfit.report_fit(res)
        self.pfit = res.params

    def mcmc(self):
        mini = lmfit.Minimizer(self._loglike, self.pfit)
        res_emcee = mini.emcee(
            burn=300, steps=600, thin=1, nwalkers=30, params=self.pfit, seed=1
        )
        self.flatchain = res_emcee.flatchain 

class FitterExponential(Fitter):
    def __init__(self, *args, **kwargs):
        super(FitterExponential, self).__init__(*args, **kwargs)
        self.model = cksmet.model.smet_exponential
        p = Parameters()
        p.add('kp',value=0.06,vary=True,min=0,max=1)
        p.add('beta',value=0.28,vary=True)
        self.p0 = p

class FitterPowerLawCutoff(Fitter):
    def __init__(self, *args, **kwargs):
        super(FitterPowerLawCutoff, self).__init__(*args, **kwargs)
        self.model = cksmet.model.powerlaw_cutoff
        p = Parameters()
        p.add('kp',value=0.06,vary=True,min=0,max=1)
        p.add('beta',value=0.28,vary=True)
        p.add('per0',value=7,vary=True,min=0,max=100)
        p.add('gamma',value=2,vary=True)
        self.p0 = p
    
