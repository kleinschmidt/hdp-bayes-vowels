from scipy.special import gammaln as lgamma
from numpy import *
import numpy as np
import numbers

class Phon:
    def __init__(self, parent=None):
        self.obs = RunningVar()
        self.count = 0
        if parent==None:
            self.params = {'nu': 1.001, 'mu': 0.0, 's': 1.0}
        else:
            self.params = parent.params
    
    def add(self, seg):
        if isinstance(seg, RunningVar):
            self.obs.merge(seg)
        elif isinstance(seg, numbers.Real):
            self.obs.push(seg)
        else:
            raise TypeError('seg needs to be number or RunningVar to add to Phon')
        self.count += 1
    
    def remove(self, seg):
        if isinstance(seg, RunningVar):
            self.obs.split(seg)
        elif isinstance(seg, numbers.Real):
            self.obs.pull(seg)
        else:
            raise TypeError('seg needs to be number or RunningVar to remove from Phon')
        self.count -= 1
    
    def lhood(self, seg, LOG=True):
        n = seg.n
        nu = self._nu_n()
        mu = self._mu_n()
        s_nu = self._s_n()
        print n, nu, mu, s_nu
        #print n, nu, mu, s_nu, lgamma, pi, log
        L = lgamma((nu+n)/2) - lgamma(nu/2) - n*0.5*log(pi*s_nu) - 0.5*log((nu+n)/nu)
        L -= (nu+n)*0.5 * log(1 + (seg.n*seg.s + (seg.m-mu)**2 * (n*nu)/(n+nu)) / s_nu)
        if not LOG:
            return exp(L)
        else:
            return L
    
    def _nu_n(self):
        return self.params['nu'] + self.obs.n
    
    def _mu_n(self):
        return (self.params['mu']*self.params['nu'] + self.obs.m*self.obs.n) / (self.params['nu']+self.obs.n)
    
    def _s_n(self):
        return (self.params['s'] + self.obs.n*self.obs.s +
                self.params['nu']*self.obs.n/(self.params['nu']+self.obs.n) * \
                    (self.obs.m - self.params['mu'])**2)

    def test(self):
        # this checks out with the R implementation in lh2.r
        p = Phon(None)
        x = range(20)
        for o in x: 
            p.add(o)
        print p.obs
        print [p.lhood(RunningVar(n=1, m=float(o), s=0.0)) for o in x]
        return p


class PhonDP:
    def __init__(self, params=None):
        if params==None:
            self.params = { 'nu': 1.001, 'mu': 0.0, 's': 1.0, 'alpha': 1.0 }
        else:
            self.params = params
        self.lexs = None
        self.phons = [Phon(self)]

    def newphon(self):
        p = Phon(self)
        self.phons.append(p)
        return p

    def iterate(self, n=1):
        for i in range(n):
            pass
        return

    def sample(self, n=1):
        probs = [p.count for p in phons]
        probs = [i/sum(probs) for i in probs]
        # multinomial returns a trials-by-phon matrix which has one
        # nonzero entry per row.  the second element returned by
        # nonzero() is the column (phon) indices.
        samp = np.nonzero(np.random.multinomial(1, probs, n))[1]
        return [self.phons[i] for i in samp]


class Lex:
    def __init__(self, parents):
        self.count = 0
        self.segs = [[RunningVar(), p] for p in parents]
    
    def add(self, word):
        for lseg, wseg in zip(self.segs, word):
            # labs[1] = parent phon
            lseg[1].add(wseg)
            # labs[0] = observations (RunningVar)
            lseg[0].push(wseg)
        self.count += 1
    
    def remove(self, word):
        for lseg, wseg in zip(self.segs, word):
            # labs[1] = parent phon
            lseg[1].remove(wseg)
            # labs[0] = observations (RunningVar)
            lseg[0].pull(wseg)
        self.count -= 1

    def holdout(self, word):
        self.remove(word)
    
    def lhood(self, word):
        L = 0
        for lseg, wseg in zip(self.segs, word):
            # lseg[1] = parent
            L += lseg[1].lhood(wseg)
        return L


class LexDP:
    def __init__(self, words):


class RunningVar:
    def __init__(self, n=0, m=0.0, s=0.0):
        self.n = n
        self.m = float(m)
        self.s = float(s)

    def __str__(self):
        return "rv(n: %d, mean: %.3f, var: %.3f)" % (self.n, self.m, self.s)
    
    def push(self, x):
        """add a single observation to this rv"""
        self.n += 1
        m_new = self.m + (x-self.m)/self.n
        self.s = ((self.n-1)*self.s + (x-self.m)*(x-m_new)) / self.n
        self.m = m_new
    
    def pull(self, x):
        """remove a single observation from this rv"""
        self.n -= 1
        if (self.n==0):
            self.m = 0
            self.s = 0
        else:
            m_new = (self.m*(self.n+1) - x) / self.n
            self.s = (self.s*(self.n+1) - (x-self.m)*(x-m_new)) / self.n
            self.m = m_new
    
    def merge(self, other):
        """merge the observations from another rv into this one"""
        if not isinstance(other, RunningVar):
            raise TypeError('merge() requires a RunningVar object')
        
        nx = self.n
        mx = self.m
        sx = self.s

        ny = other.n
        my = other.m
        sy = other.s

        self.n = nx + ny
        self.m = (nx*mx + ny*my)/self.n
        self.s = (nx*sx + ny*sy)/self.n + (mx-my)**2 * nx*ny/(nx+ny)**2
    
    def split(self, other):
        """remove multiple observations from this rv"""
        if not isinstance(other, RunningVar):
            raise TypeError('split() requires a RunningVar object')
        
        n = self.n
        m = self.m
        s = self.s
        
        self.n -= other.n
        self.m = (n*m - other.n*other.m) / self.n
        self.s = s*n/self.n - (self.m-other.m)**2 * other.n/n - other.n*other.s / self.n
    
    def rvTest(self):
        """test each of the four runningVar methods"""
        rv1 = RunningVar()
        rv2 = RunningVar()
        x1 = range(0,10)
        x2 = range(9,20)
        for x in x1:
            rv1.push(x)
        for x in x2:
            rv2.push(x)
        print 'test of rv.push()'
        print rv1
        print meanVar(x1)
        print 'test of rv.pull()'
        rv1.pull(9)
        print rv1
        print meanVar(range(0,9))
        print 'test of rv.merge()'
        rv1.merge(rv2)
        print rv1
        print meanVar(range(0,20))
        print 'test of rv.split()'
        rv3 = RunningVar()
        for x in x1:
            rv3.push(x)
        rv1.split(rv3)
        print rv1
        print meanVar(range(10,20))
    
    def meanVar(xx):
        """offline mean/variance computation to compare with online"""
        n = len(xx)*1.0
        m = sum(xx)/n
        s = sum([(x-m)**2 for x in xx])/n
        return m,s


################################ TESTING STUFF ##########################################
def makeGaussianPhon(mean=0.0, var=1.0, n=20):
    samp = np.random.normal(mean, var, n)
    p = Phon()
    for x in samp: p.add(x)
    return p, samp

def makeSomeWords():
    obs1 = np.random.normal(2, 1.0, 50)
    obs2 = np.random.normal(-2, 1.0, 50)
    words1 = zip(obs1[0:25], obs2[0:25])
    words2 = zip(obs2[25:50], obs1[25:50])
    return (words1, words2)

