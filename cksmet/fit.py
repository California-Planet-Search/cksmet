import numpy as np
from lmfit import Parameters 
import lmfit
import cksmet.model
import copy
import corner
from lmfit import Parameters
import numpy as np

fmt = {'kp':"{:.3f}",'logkp':"{:.2f}",'beta':"{:+.1f}",'alpha':"{:+.1f}",'gamma':"{:.1f}",'per0':"{:.1f}",'binwper':"{:.1f}"}

class Fit(object):
    def __init__(self, x, nplnt, ntrial):
        """
        x (array): centers of bins (can be 1D or 2D)
        nplnt (array): number of planets per bin
        ntrial (array): number of trials per bin
        """
        
        self.x = np.array(x)
        self.nplnt = np.array(nplnt)
        self.ntrial = np.array(ntrial)

    def _loglike(self, params):
        """
        Args:
            params (lmfit.Parameters): parameters
        """
        F = self.model(params, self.x) # Predicted occurrence at each bin.
        F = F.reshape(*self.nplnt.shape)
        if (F > 1).any():
            return -1e10 # doesn't work for negative likelihood

        nnd = self.ntrial - self.nplnt # number of non detections (bin by bin)
        _loglike = np.sum(self.nplnt * np.log(F)) + np.sum(nnd * np.log(1 - F) )
        if np.isnan(_loglike):
            return -1e10

        return _loglike

    def _negloglike(self,params):
        return -1.0 * self._loglike(params)

    def print_parameters(self,prefix=''):
        lines = self.to_string()
        print "\n".join(lines)
        
    def to_string(self,prefix=''):
        lines = []
        for k in self.pfit.keys():
            s = fmt[k]
            if self.pfit[k].vary:
                chain = self.flatchain[k]
                q = chain.quantile([0.16,0.50,0.84])
                val = s.format(q.ix[0.50])
                err1 = s.format(q.ix[0.84] - q.ix[0.50])
                err2 = s.format(q.ix[0.16] - q.ix[0.50])
                s = "$%s^{+%s}_{%s}$" % (val, err1,err2)
                s = s.replace('++','+')
            else:
                s = "$%s$" % s.format(self.pfit[k].value)

            line = r"{{{}{}}}{{{}}}".format(prefix,k,s)
            lines.append(line)

        return lines
        
    def sample_chain(self, nsamp):
        psamples = self.flatchain.sample(nsamp)
        for k in self.p0.valuesdict().keys():
            if not self.p0[k].vary:
                psamples[k] = self.p0[k]
        return psamples

    def sample(self, x, nsamp,):
        # method only works for 1-D models        
        psamples = self.sample_chain(nsamp)
        
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
        self.loglike_fit = self._loglike(self.pfit)

    def mcmc(self, **kwargs):
        mini = lmfit.Minimizer(self._loglike, self.pfit, nan_policy='raise')
        res_emcee = mini.emcee(params=self.pfit, seed=1, **kwargs)
        self.flatchain = res_emcee.flatchain 
        
        print "loglike_fit = {:.4f}, max(loglike_mcmc) = {:.4f}".format(self.loglike_fit, np.max(res_emcee.lnprob) )
        if self.loglike_fit < np.max(res_emcee.lnprob):
            print "warning maxlike from lmfit is smaller than maxlike from emcee"

    def corner(self):
        corner.corner(self.flatchain)

class Exponential(Fit):
    def __init__(self, *args, **kwargs):
        super(Exponential, self).__init__(*args, **kwargs)
        self.model = smet_exp
        p = Parameters()
        p.add('logkp',value=-2,vary=True,min=-5,max=1)
        p.add('beta',value=0.28,vary=True,min=-10,max=10)
        self.p0 = p

class PerPowerLawExpSmetExp(Fit):
    def __init__(self, *args, **kwargs):
        super(PerPowerLawExpSmetExp, self).__init__(*args, **kwargs)
        self.model = per_powerlaw_smet_exp
        p = Parameters()
        p.add('logkp',value=-2,vary=True,min=-5,max=1)
        p.add('alpha',value=1.8,vary=False,min=-6,max=6)
        p.add('beta',value=0.28,vary=True,min=-6,max=6)
        p.add('binwper',vary=False,min=-6,max=6)
        self.p0 = p


class PowerLawCutoff(Fit):
    def __init__(self, *args, **kwargs):
        super(PowerLawCutoff, self).__init__(*args, **kwargs)
        self.model = powerlaw_cutoff
        p = Parameters()
        p.add('logkp',value=-2,vary=True,min=-5,max=1)
        p.add('beta',value=0.28,vary=True)
        p.add('per0',value=7,vary=True,min=0,max=100)
        p.add('gamma',value=2,vary=True)
        self.p0 = p

def powerlaw_cutoff(params, x):
    """
    Equation 8 in Howard 2010
    """
    p = params.valuesdict()
    per = x
    occur = (
        10**p['logkp']  * 
        per**p['beta'] * 
        (1 - np.exp( -1.0 * (per / p['per0'])**p['gamma']))
     ) 
    return occur

def smet_exp(params, x):
    """
    df/dM = kp * 10**(beta * M)
    """
    smet = x
    p = params.valuesdict()
    occur = 10**p['logkp'] * 10**(p['beta'] * smet)
    return occur
    
def per_powerlaw_smet_exp(params, x):
    """
    df / dM dlogP = kp * P**alpha 10**(beta * M)
    """
    per = x[:,0]
    smet = x[:,1]
    p = params.valuesdict()
    
    # occur is number of planets per star per log period interval
    occur = 10**p['logkp'] * per**p['alpha'] * 10**(p['beta'] * smet)

    # Convert into number of planets per bin (compared to data)
    occur *= p['binwper']
    return occur


from scipy.special import xlogy
def per_powerlaw_smet_exp_integral(p, lims):
    """
    Evaluate double integral over [per1, per2] and [smet1, smet2]
    """
    per1 = lims[0][0]
    per2 = lims[0][1]
    smet1 = lims[1][0]
    smet2 = lims[1][1]
    
    # Integral is per**alpha / alpha. This gives problems for alpha=0
    int1 = lambda per: per** p['alpha'] / p['alpha']
    int2 = lambda smet : 10**(p['beta'] * smet) / p['beta'] / np.log(10)
    if p['alpha']==0:
        f = 10**p['logkp']*(int2(smet2) - int2(smet1))
    else:
        f = 10**p['logkp']*(int1(per2) - int1(per1))*(int2(smet2) - int2(smet1))

    return f

