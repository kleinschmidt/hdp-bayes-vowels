import numpy as np
from scipy.special import multigammaln as lgammad

def onlineVar(data):
    M2 = 0.
    n = 0
    mean = 0.

    for x in data:
        n += 1
        delta = x-mean
        mean += delta/n
        M2 += np.outer(delta, (x-mean))

    return n, mean, M2/n


def mvlhood(obs, cat, LOG=True):
    nu0 = 1.001
    mu0 = zeros(cat.mean.shape)
    s0 = diag(ones(cat.mean.shape), 0)
    
    nu = nu0 + cat.n
    mu = (mu0*nu0 + cat.n*cat.mean) / nu
    s = s0 + cat.M2 + np.outer(cat.mean-mu0, cat.mean-mu0) * nu0*cat.n/nu
    d = max(cat.mean.shape)
    print nu, mu, s, d

    lgam = lambda x: lgammad(x, d)
    L = lgam((nu+obs.n)/2) + nu/2 * np.log(np.linalg.det(s)) - \
        lgam(nu/2) - d*obs.n/2 * np.log(np.pi) - d/2 * np.log(obs.n/nu + 1) - \
        (nu+obs.n)/2 * np.linalg.det(s + obs.M2 + np.outer(obs.mean-mu, obs.mean-mu) * obs.n*nu/(obs.n+nu))
    if LOG: return L
    else: return np.exp(L)

def mvlhood1(obs, cat, LOG=True):
    return mvlhood(MDRunningVar(n=1, mean=obs, M2=0.), cat, LOG)

class MDRunningVar():
    def __init__(self, n=0, mean=0., M2=0.):
        self.n = n          # count of observations
        self.mean = mean    # mean of observations
        self.M2 = M2        # second cumulant (\sum_i (x_i-\bar{x})^2)

    def __str__(self):
        return "mrv(n=%s, m=%s, s=%s)" % (self.n, self.mean, self.var())

    def __repr__(self):
        return "mrv(n=%s, m=%s, s=%s)" % (repr(self.n), repr(self.mean), repr(self.var()))
        
    def var(self):
        """Convert second cumulant (which is tracked incrementally) into variance"""
        if self.n > 0: return self.M2 / self.n
        else: return 0

    def push(self, x):
        delta = x - self.mean
        self.n += 1
        self.mean += delta/self.n
        self.M2 += np.outer(delta, (x-self.mean))

    def pushall(self, X):
        for x in X: self.push(x)

    def pull(self, x):
        delta = x - self.mean
        self.n -= 1
        self.mean -= delta/self.n
        self.M2 -= np.outer((x-self.mean), delta)

    def pullall(self, X):
        for x in X: self.pull(x)

    def merge(self, other):
        if not isinstance(other, MDRunningVar):
            raise TypeError('Can only merge another MDRunningVar')
        delta = other.mean - self.mean
        M2delta = np.outer(delta, delta) * self.n*other.n / (self.n+other.n)
        self.n += other.n
        self.mean += delta*other.n/self.n
        self.M2 += other.M2 + M2delta

    def split(self, other):
        if not isinstance(other, MDRunningVar):
            raise TypeError('Can only split off another MDRunningVar')
        self.n -= other.n
        self.mean -= (other.mean-self.mean) * other.n/self.n
        self.M2 -= other.M2 + np.outer(self.mean-other.mean, self.mean-other.mean) * self.n*other.n/(self.n+other.n)
