from lmfit import Parameters
import numpy as np


# Multiple broken powerlaw
'''
# Must be in order
p = Parameters()
p.add('per1',value=1,vary=False)
p.add('per2',value=10,vary=True)
p.add('per3',value=300,vary=False)
p.add('logf1',value=0,vary=True)
p.add('logf2',value=1,vary=True)
p.add('logf3',value=1,vary=True)

def occur(p, peri):
    """
    Return log10 of occurrence per bin 
    """
    pdict = p.valuesdict()
    per = [p[k].value for k in pdict.keys() if k.count('per')==1]
    per = np.array(per)
    logper = np.log10(per)
    logf = [p[k].value for k in pdict.keys() if k.count('logf')==1]
    logperi = np.log10(peri)
    return np.interp(logperi,logper,logf)

clf()
peri = logspace(log10(1),log10(00),100)
semilogx()
plot(x,occur(p, x ))

cut = df[(df.prad2==1.7) & (df.smet2==0.0) & (df.prob_det_mean > 0.1)]
perc = cut.perc
f = cut.nplnt / cut.ntrial

def obj(p, perc, f):
    """
    Chisq is the difference (logf_model - logfobs)**2
    """
    logf_model = occur(p,perc)
    logf_data = np.log10(f)
    _resid = logf_model - logf_data 
    return _resid

res = lmfit.minimize(obj, p, args=(perc,f))
print fit_report(res)
clf()
loglog()
plot(cut.perc,cut.nplnt / cut.ntrial,'.')
plot(peri,10**occur(res.params , peri))

'''


# Knee and slope

def per_occur_knee_slope(params, per):
    """
    Toy model for occurrence which is define by a knee and two slopes
    coming off it

    logf = 
         m1 (logper - logperm) + logfm # if logper < logperm
         m2 (logper - logperm) + logfm # if logper > logperm

    """
    
    logper = np.log10(per)
    p = params.valuesdict()

    b1 = logper < p['logperm']
    b2 = logper > p['logperm']
    logf = np.empty_like(per)
    
    # Left side
    logf[b1] = p['m1']  * (logper[b1] - p['logperm']) + p['logfm']
    logf[b2] = p['m2']  * (logper[b2] - p['logperm']) + p['logfm']
    return logf


def powerlaw_and_cutoff(params, per):
    """
    Equation 8 in Howard 2010
    """
    p = params.valuesdict()
    occur = (
        p['kp']  * 
        per**p['beta'] * 
        (1 - np.exp( -1.0 * (per / p['per0'])**p['gamma']))
     ) 
    return occur
    
def loglike_powerlaw_and_cutoff(params, perc, nplnt, ntrial):
    F = powerlaw_and_cutoff(params, perc) # Predicted occurrence at each bin.
    nnd = ntrial - nplnt # number of non detections (bin by bin)
    _loglike = np.sum(nplnt * np.log(F)) + np.sum(nnd * np.log(1 - F) )
    return _loglike
