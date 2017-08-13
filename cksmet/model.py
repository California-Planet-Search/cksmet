from lmfit import Parameters
import numpy as np

def powerlaw_cutoff(params, per):
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

def smet_exponential(params, smet):
    """
    Equation 8 in Howard 2010

    df / d [Fe/H] = kp * 10**(beta * [Fe/H])
    """
    p = params.valuesdict()
    occur = p['kp'] * 10**(p['beta'] * smet)
    return occur
    
def loglike_1d(params, func, x, nplnt, ntrial):
    """
    Args:
        params (lmfit.Parameters): parameters
        func (function): occurrence model. First arg must be params, second 
            must be independent variable
        x (array): independent var
        nplnt (array): number of detections
        ntrial (array): number of non-detections
    """
    F = func(params, x) # Predicted occurrence at each bin.
    nnd = ntrial - nplnt # number of non detections (bin by bin)
    _loglike = np.sum(nplnt * np.log(F)) + np.sum(nnd * np.log(1 - F) )
    return _loglike

