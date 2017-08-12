"""
Module for statistics
"""
class Binomial(object):
    def __init__(self, n, k):
        """
        Args:
            n: number of trials
            k: number of sucesses
        """
        self.n = float(n)
        self.k = float(k)
        self.npdf = 10000
        if k!=0:
            # Detection, calculate pdf of the rate of planets in box
            rate = 1.0 * k / n
            self.maxrate = min(1,10*rate)
        else:
            # Non-detection, upper limit on rate.
            self.maxrate = min(1, 10 / self.n)

    def pdf(self):
        """
        Returns probabilty of rate r
        """
        _rate = np.linspace(0,self.maxrate,self.npdf)
        _pdf = binom_gamma(_rate,self.n,self.k)
        _pdf /= _pdf.sum()
        return _rate, _pdf

    def sample(self, nsamp):
        rate, pdf = self.pdf()
        cdf = pdf.cumsum()
        fp = rate
        x = np.random.random_sample(nsamp)
        return interp(x, cdf, rate)

    def hist_samples(self, nsamp, downsamp=10):
        """
        Verify that the sampling working
        """
        
        rate, pdf = self.pdf()
        samples = self.sample(nsamp)
        weights = ones_like(samples)/nsamp # normalize histogram
        hist(samples,bins=rate[::downsamp],weights=weights)
        plot(rate[::downsamp],pdf[::downsamp]*downsamp)

def binom_gamma(p, n, k):
    """
    Probability that rate is p given there are n trials and k sucesses
    """
    n = 1.0*n
    k = 1.0*k
    p = 1.0*p
    combiln = (gamln(n+1) - (gamln(k+1) + gamln(n - k + 1)))
    _logpdf = combiln + special.xlogy(k, p) + special.xlog1py(n-k, -p)
    _pdf = np.exp(_logpdf)
    return _pdf

def sum_cells(ntrial, nplnt):
    """
    Add the occurrence values from different cells.

    Args:
        nplnt (arary): number of planets detected
    """

    ntrial = np.array(ntrial)
    nplnt = np.array(nplnt)
    nsample = int(1e4)

    samplesL = []
    for i in range(len(nplnt)):
        binom = Binomial(ntrial[i], nplnt[i])
        samplesL.append(binom.sample(nsample))

    samplesL = np.vstack(samplesL)
    isuplim = (nplnt==0) # True if cell yields upper limit
    

    samples_sum = samplesL[~isuplim].sum(0)
    d = {}
    d = dict(rate=None, rate_err1=None, rate_err2=None, rate_ul=None)

    # All measurements are upper limits
    if (nplnt==0).all():
        samples_sum = samplesL.sum(0)
        p16, p50, p84, p90 = np.percentile(samples_sum, [16,50,84,90])
        d['rate_ul'] = p90
    else:
        samples_sum = samplesL[~isuplim].sum(0)
        p16, p50, p84, p90 = np.percentile(samples_sum, [16,50,84,90])
        d['rate'] = p50
        d['rate_err1'] = p84 - p50
        d['rate_err2'] = p16 - p50
    return d

