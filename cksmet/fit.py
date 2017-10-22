import numpy as np
from lmfit import Parameters 
import lmfit
import cksmet.model
import copy
import corner
from lmfit import Parameters
import numpy as np

fmt = {'kp':"{:.3f}",'logkp':"{:.2f}",'beta':"{:+.1f}",'alpha':"{:+.1f}",'gamma':"{:.1f}",'gamma+beta':"{:.1f}",'per0':"{:.1f}"}

class Fit(object):
    def __init__(self, x, dx, nplnt, ntrial):
        """
        x (array): centers of bins (can be 1D or 2D)
        dx (array): volume of box
        nplnt (array): number of planets per bin
        ntrial (array): number of trials per bin
        """
        self.x = np.array(x)
        self.dx = dx
        self.nplnt = np.array(nplnt)
        self.ntrial = np.array(ntrial)

    def _loglike(self, params):
        """
        Args:
            params (lmfit.Parameters): parameters
        """
        # Planet distribution function evaluated at the center of each cell
        dfdx = self.model(params, self.x) 
        
        # Number of expected planet detections from model integrated
        # over each cell
        fcell = dfdx * self.dx
        fcell = fcell.reshape(*self.nplnt.shape)
        if (fcell > 1).any():
            print "fcell > 1" 
            return -1e10 # doesn't work for negative likelihood

        nnd = self.ntrial - self.nplnt # number of non detections (bin by bin)
        _loglike = np.sum(
            self.nplnt * np.log(fcell)) + np.sum(nnd * np.log(1 - fcell) 
        )
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
        p.add('logkp',value=-2,vary=True,min=-10,max=10)
        p.add('beta',value=0.28,vary=True,min=-10,max=10)
        self.p0 = p

class PerPowerLawExpSmetExp(Fit):
    def __init__(self, *args, **kwargs):
        super(PerPowerLawExpSmetExp, self).__init__(*args, **kwargs)
        self.model = per_powerlaw_smet_exp
        p = Parameters()
        p.add('logkp',value=-2,vary=True,min=-10,max=10)
        p.add('alpha',value=1.8,vary=False,min=-6,max=6)
        p.add('beta',value=0.28,vary=True,min=-6,max=6)
        self.p0 = p

class PowerLawCutoff(Fit):
    def __init__(self, *args, **kwargs):
        super(PowerLawCutoff, self).__init__(*args, **kwargs)
        self.model = powerlaw_cutoff
        p = Parameters()
        p.add('logkp',value=0,vary=True,min=-10,max=10)
        p.add('beta',value=0.28,vary=True)
        p.add('per0',value=7,vary=True,min=0,max=100)
        p.add('gamma',value=2,vary=True)
        self.p0 = p

    def to_string(self,prefix=''):
        # Prints the sum of gamma+beta
        lines = super(PowerLawCutoff,self).to_string(prefix=prefix)
        k =  'gamma+beta'
        chain = self.flatchain['gamma'] + self.flatchain['beta'] 
        s = fmt[k]
        q = chain.quantile([0.16,0.50,0.84])
        val = s.format(q.ix[0.50])
        err1 = s.format(q.ix[0.84] - q.ix[0.50])
        err2 = s.format(q.ix[0.16] - q.ix[0.50])
        s = "$%s^{+%s}_{%s}$" % (val, err1,err2)
        s = s.replace('++','+')
        line = r"{{{}{}}}{{{}}}".format(prefix,k,s)
        lines.append(line)
        return lines

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
    return occur


from scipy.special import xlogy

def per_powerlaw_smet_exp_integral(p, lims, smet):
    """
    Equation 15a from Youdin 11

    Args:
        p: parmaeter object
        lims: period integration limits
        smet: metalicity of stars

    For each star, 
    1. evaluate model at df/dx(per,smet),
    2. integrate over orbital period
    3. 
    
    """
    per1 = lims[0]
    per2 = lims[1]
    nstars = len(smet)
    def _int(per):
        if p['alpha']==0:
            return np.log(per) / np.log(10)
        else:
            return (per**p['alpha']) / p['alpha'] / np.log(10)

    kp = 10**p['logkp']
    f = kp / nstars * np.sum( 10**(p['beta']*smet) ) * (_int(per2) - _int(per1))
    return f

