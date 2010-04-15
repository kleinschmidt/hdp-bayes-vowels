from scipy.special import gammaln as lgamma
from numpy import *
import numpy as np
import numbers, warnings, sys


def debugMethodInfo(func):
    def wrapper(self, *args, **kwargs):
        sys.stdout.write( '%s.%s ' % (self.__class__.__name__, func.__name__) )
        func(self, *args, **kwargs)
        sys.stdout.write( ' (done)\n' )
    return wrapper


class PhonDP:
    def __init__(self, params=None):
        if not params:
            self.params = { 'nu': 1.001, 'mu': 0.0, 's': 1.0, 'alpha': 1.0 }
        else:
            self.params = params
        self.children = None
        self.phoncounter = 0
        self.phons = []
        self.newphon()
        #self.priors = [Phon(parent=self, count=self.params['alpha']/self.params['s']) for n in range(self.params['s'])]
        # in a situation where the prior wasn't simple the above expression should
        # be used (so the prior probability is approximated by sampling
        # s items from the base distribution) but since the prior in this case is
        # analytically tractable (with an empty Phon) there's no need to bother
        self.priors = [Phon(parent=self, count=self.params['alpha'])]
    
    def newphon(self):
        p = Phon(self, id=self.phoncounter)
        self.phoncounter += 1
        self.phons.append(p)
        print '  adding Phon', p.id
        return p
    
    def iterate(self, n=1):
        print 'Phon sweep...'
        for i in range(n):
            # iterate n times
            for seg in self.children.segments():
                # seg should be a Segment instance, check here
                if not isinstance(seg, Segment): raise TypeError('seg should be a Segment')
                # hold this segment out
                #try:
                lab = seg.holdout()
                #except EmptyTable:
                if lab.count == 0:
                    #...pruning its label if necessary
                    self.prune(lab)
                # calculate probabilities (lhood + CRP prior) and sample a Phon
                probs = [np.exp(p.lhood(seg.obs)) * p.count for p in self.phons+self.priors]
                i = sampleMultinomial(probs)
                # record the new Phon, creating a new one if necessary
                if i < len(self.phons):
                    # picked existing phon
                    seg.relabel(self.phons[i])
                else:
                    # picked new phon
                    seg.relabel(self.newphon())
        return

    def prune(self, phon):
        """
        Remove a given phon from this DP (method exists pretty much only for
        symmetry's sake with LexDP
        """
        print '  pruning Phon', phon.id
        self.phons.remove(phon)
    
    def sample(self, n=1):
        """
        Sample from this DP.  Returns a list of n Phons, each with probability 
        weighted by their count in the DP
        """
        probs = [p.count for p in self.phons]
        return [self.phons[i] for i in sampleMultinomial(probs, n)]


class LexDP:
    def __init__(self, words, parent=None, params=None):
        # initialize params
        if not params:
            self.params = { 'nu': 1.001, 'mu': 0.0, 's': 1.0, 'beta': 1.0, 'alpha': 1.0, 'r': 5 }
        else:
            self.params = params
        # initialize parent PhonDP
        if not parent:
            self.parent = PhonDP(params)
        else:
            self.parent = parent
        self.parent.children = self
        self.lexcounter = 0
        # sort words by length
        wordsByLen = dict()
        for w in words:
            try:
                wordsByLen[len(w)] += [w]
            except KeyError:
                wordsByLen[len(w)] = [w]
        # initialize lexs with one Lex per unique word length
        self.lexs = []
        self.priors = []
        self.words = []
        for length in wordsByLen.keys():
            # create and add a new Lex for this length
            newlex = self.addLex(length=length)
            # add each word of the right length to it (and it's Phons...) and
            # append a properly labeled Word to self.words
            for w in wordsByLen[length]:
                newlex.add(w)
                self.words.append(Word(lex=newlex, obs=w))
            # initialize priors, too
            self.priors.extend(self.sampleBase(length=length, 
                                               n=self.params['r'], 
                                               count=self.params['beta']/self.params['r']))

    def __str__(self):
        retr = 'LexDP state:\n' + \
            ' Labels:\n  ' + str([w.lex.id for w in self.words]) + '\n' + \
            ' Lexicon:\n'
        for lex in self.lexs:
            retr += '  %s -> (n=%d) %s\n' % (lex.id, lex.count, [s.phon.id for s in lex])
        #retr += ' Priors:\n'
        #for lex in self.priors:
        #    retr += '  %s -> (n=%d) %s\n' % (lex.id, lex.count, [s.phon.id for s in lex])
        retr += ' Phones:'
        for p in self.parent.phons:
            retr += '\n  %s -> (n=%d) %s' % (p.id, p.count, p.obs)
        return retr
    
    def iterate(self, n=1):
        for i in range(n):
            print 'Lex sweep...'
            for word in self.words:
                # hold out this word
                oldlex = word.holdout()
                # pruning it's Lex label if necessary
                if oldlex.count == 0: self.prune(oldlex)
                # calculate probability of each lex (log lhood + log (pseudo)count)
                probs = [np.exp(lex.lhood(word)) * lex.count for lex in self.lexs+self.priors]
                i = sampleMultinomial(probs)
                # record the new label, adding a new lex/refreshing priors as necessary
                if i < len(self.lexs):
                    # old lex
                    word.relabel(self.lexs[i])
                else:
                    # new lex...get sampled one from prior and replace it...
                    newlex = self.priors[i-len(self.lexs)]
                    # ...and add it to the list of Lexs
                    self.addLex(newlex)
                    # record the new Lex
                    word.relabel(newlex)
                    # refresh priors (since phon segment counts have changed)
                    self.refreshPriors()
            # sweep through the parent PhonDP
            self.parent.iterate()
            self.refreshPriors()
            # output the current state just to see what's up
            print self

    def refreshPriors(self):
        """
        Resample Lexs from the base distribution to be used in approximation
        of the prior probability of a word.  This should be done whenever the
        counts of any of the Phons in the parent distribution change.
        """
        self.priors = [self.sampleBase(length=len(lex),
                                       count=self.params['beta']/self.params['r']).pop()
                       for lex in self.priors]
    
    def __new_OLD__refreshPriors(self, length=None):
        """
        Resample Lexs of the specified length from the base distribution defined
        the parent PhonDP.
        """
        # get the indices to replace--either the indices w/ len==length or all
        # the priors if no length is specified
        indices = (i for i in range(len(self.priors))
                   if len(self.priors[i])==length or not length)
        print indices
        for i in indices:
            self.priors[i] = self.sampleBase(len(self.priors[i]), n=1,
                                             count=self.params['beta']/self.params['r']).pop()
    
    def __OLD__refreshPriors(self, length=None):
        for l in self.priors:
            if not length or len(l)==length:
                self.priors.remove(l)
                self.priors.extend(self.sampleBase(len(l), n=1,
                                                   count=self.params['beta']/self.params['r']))
    
    def prune(self, lex):
        """Remove a given Lex from the list, pruning parent Phons as necessary"""
        print '  pruning Lex %s (n=%d)' % (lex.id, lex.count)
        for seg in lex:
            # throw in a sanity check, just for the hell of it.  if lex.count==0
            # then seg.count should == 0 for all segments.
            if seg.obs.n != lex.count:
                raise Exception("Segment (%d) and Lex (%d) counts disagree..." % 
                                (seg.obs.n, lex.count))
            # remove from parents, too
            #try:
            #    seg.phon.remove(seg)
            #except EmptyTable:
            #    self.parent.prune(seg.phon)
            seg.phon.remove(seg)
            if seg.phon.count == 0:
                self.parent.prune(seg.phon)
        self.lexs.remove(lex)
        self.refreshPriors()
    
    def addLex(self, newlex=None, length=None):
        """Add a new Lex to this DP, sampling the base if necessary"""
        if not newlex:
            if not length:
                raise ValueError('Must specify a length for Lex to be sampled')
            newlex = self.sampleBase(length, n=1, count=0).pop()
        for seg in newlex: seg.phon.add(seg)
        # initialize lex.count as well (only necessary if newlex was a prior before...)
        newlex.count = 0
        newlex.id = self.lexcounter
        self.lexcounter += 1
        self.lexs.append(newlex)
        # print some output for debugging
        print '  adding Lex %s: %s' % (newlex.id, [seg.phon.id for seg in newlex])
        return newlex
    
    def sampleBase(self, length, n=1, count=0):
        """Sample n Lexs from the base distribution defined by self.parent"""
        return [Lex(self.parent.sample(length), count) for i in range(n)]
    
    def sample(self, n=1, length=None):
        """
        Sample a Lex from this DP, with probability proportional to each Lex's
        count.

        Keyword Args:
        n = 1 (number of samples)
        length = None (if non-none, only samples Lex's of given length)
        """
        if not length:
            probs = [l.count*(len(l)==length) for l in self.lexs]
        else:
            probs = [l.count for l in self.lexs]
        return [self.lexs[i] for i in sampleMultinomial(probs, n)]
    
    def segments(self):
        """Generates an iterator over all segments in the lexicon"""
        for lex in self.lexs:
            for seg in lex:
                yield seg


class Phon:
    def __init__(self, parent=None, count=0, id=None):
        self.obs = RunningVar()
        self.count = count
        self.id = id
        if parent==None:
            self.params = {'nu': 1.001, 'mu': 0.0, 's': 1.0}
        else:
            self.params = parent.params

    def __str__(self):
        return "Phon %s: count=%d, obs=%s" % (self.id, self.count, self.obs)

    def add(self, seg):
        """
        Add a Segment to this Phon, updating summary stats and .count
        """
        if not isinstance(seg, Segment):
            raise TypeError('seg needs to be a Segment instance to add to Phon (not %s)' % seg.__class__)
        self.obs.merge(seg.obs)
        self.count += 1
    
    def remove(self, seg):
        """
        Remove a segment from this Phon, updating summary stats and .count
        """
        if not isinstance(seg, Segment):
            raise TypeError('seg needs to be a Segment instance to remove from Phon (not %s)' % seg.__class__)
        self.obs.split(seg.obs)
        self.count -= 1
#        if self.count == 0: raise EmptyTable('Prune this empty Phon!')
        if self.count < 0: raise ValueError('Phon count is less than 0...')

    def push(self, seg):
        """
        Push additional observations into this Phon, updating summary stats but
        not .count (use .add() instead)
        """
        if isinstance(seg, RunningVar):
            self.obs.merge(seg)
        elif isinstance(seg, numbers.Real):
            self.obs.push(seg)
        else:
            raise TypeError('seg needs to be number or RunningVar to push to Phon (not %s)' % seg.__class__)

    def pull(self, seg):
        """
        Pull some observations from this Phon, updating summary stats but not
        .count (use .remove() instead)
        """
        if isinstance(seg, RunningVar):
            self.obs.split(seg)
        elif isinstance(seg, numbers.Real):
            self.obs.pull(seg)
        else:
            raise TypeError('seg needs to be number or RunningVar to pull from Phon (not %s)' % seg.__class__)        
        
    def lhood(self, seg, LOG=True):
        """
        Calcunlate the likelihood (defaults to log) of a given segment, which should
        be a RunningVar instance.
        """
        if not isinstance(seg, RunningVar):
            raise TypeError('seg must be a RunningVar instance (not %s)' % seg.__class__)
        n = seg.n
        nu = self._nu_n()
        mu = self._mu_n()
        s_nu = self._s_n()
        #print n, nu, mu, s_nu
        #print n, nu, mu, s_nu, lgamma, pi, log
        L = lgamma((nu+n)/2) - lgamma(nu/2) - n*0.5*log(pi*s_nu) - 0.5*log((nu+n)/nu)
        L -= (nu+n)*0.5 * log(1 + (seg.n*seg.s + (seg.m-mu)**2 * (n*nu)/(n+nu)) / s_nu)
        if not LOG:
            return exp(L)
        else:
            return L

    def lhood1(self, seg, LOG=True):
        """
        Calculate the likelihood (defaults to log-lhood) for a single observation,
        which is either numeric or an ndarray.

        For likelihood of a group of observations, use .lhood().
        """
        nu = self._nu_n()
        mu = self._mu_n()
        s_nu = self._s_n()
        #print nu, mu, s_nu
        #print n, nu, mu, s_nu, lgamma, pi, log
        L = lgamma((nu+1)/2) - lgamma(nu/2) - 0.5*log(pi*s_nu) - 0.5*log((nu+1)/nu)
        L -= (nu+1)*0.5 * log(1 + (seg-mu)**2 * nu/(1+nu) / s_nu)
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
        x = [Segment(RunningVar(n=1,m=o,s=0.), self) for o in range(20)]
        for o in x: 
            self.add(o)
        print self.obs
        print [self.lhood(o.obs) for o in x]
        return self


class Lex:
    def __init__(self, parents, count=0, id=None):
        self.count = count
        self.segs = [Segment(RunningVar(), p) for p in parents]
        self.id = id
    
    def __iter__(self):
        return self.segs.__iter__()

    def __str__(self):
        if self.id == None:
            ret = 'Lex: unnamed'
        else:
            ret = 'Lex: ' + str(self.id)
        for seg in self:
            ret += '\n  ' + seg.__str__()
        return ret
    
    def __len__(self):
        """Return the length of this Lex in segments"""
        return len(self.segs)
    
    def phons(self):
        """Return this Lex's Phons"""
        return [s.phon for s in self.segs]
    
    def obs(self):
        """Return the obs (RunningVar) for each of this Lex's segments"""
        return [s.obs for s in self.segs]
    
    def add(self, word):
        """
        Add a single word to this Lex (and the corresponding segments from
        the corresponding Phons)
        """
        for lseg, wseg in zip(self.segs, word):
            lseg.phon.push(wseg)
            lseg.obs.push(wseg)
        self.count += 1
    
    def remove(self, word):
        """
        Remove a single word from this Lex (and the corresponding segments
        from the corresponding Phons
        """
        for lseg, wseg in zip(self.segs, word):
            lseg.phon.pull(wseg)
            lseg.obs.pull(wseg)
        self.count -= 1
#        if self.count == 0: raise EmptyTable('Prune this empty Lex!')
        if self.count < 0: raise ValueError('Lex count should not be less than 0...')
    
#    def holdout(self, word):
#        self.remove(word)
    
    def lhood(self, word):
        if len(self) != len(word):
            return -np.inf
        else:
            L = 0
            for lseg, wseg in zip(self.segs, word.obs):
                L += lseg.phon.lhood1(wseg)
            return L


class Word:
    """
    Simple container class, to bind together a list of observations and a label.

    Fields:
    obs (should be a list/tuple of numbers or ndarrays)
    lex (should be a Lex instance)
    """
    def __init__(self, obs, lex):
        if not isinstance(lex, Lex):
            raise TypeError('Word lex label must be a Lex instance (not %s)' % lex.__class__)
        self.obs = obs
        self.lex = lex

    def __iter__(self):
        return self.obs.__iter__()

    def __len__(self):
        return len(self.obs)

    def relabel(self, newlex):
        """Apply a new Lex label to this word"""
        self.lex = newlex
        self.lex.add(self.obs)
    
    def holdout(self):
        """Remove this word from its Lex, set lex to None, and return the old Lex"""
        l = self.lex
        self.lex = None
        l.remove(self.obs)
        return l


class Segment:
    """
    Simple container class, to bind together a single observation and a label.

    Fields:
    obs (should be a RunningVar instance)
    phon (should be a Phon instance)
    """
    def __init__(self, obs, phon):
        if not isinstance(obs, RunningVar):
            raise TypeError('Segment obs should be a RunningVar instance')
        if not isinstance(phon, Phon):
            raise TypeError('Segment phon must be a Phon instance')
        self.obs = obs
        self.phon = phon

    def __str__(self):
        return "Segment: %s, %s" % (self.phon.__str__(), self.obs.__str__())
    
    def relabel(self, newphon):
        """Apply a new label to this segment"""
        self.phon = newphon
        newphon.add(self)
    
    def holdout(self):
        """Hold this segment out of its Phon, set phon to None, and return the old label"""
        self.phon.remove(self)
        p = self.phon
        self.phon = None
        return p


class EmptyTable(Exception):
    """
    A custom Exception to indicate that the last customer has been removed from
    a table (Lex of Phon).  This should be used to trigger some sort of pruning
    further up the stack (in something like .holdout() or .iterate())
    """
    def __init__(self, value):
        self.value = value
        print str(value)
    def __str__(self):
        return repr(self.value)


class RunningVar:
    """
    RunningVar tracks the mean, variance, and count of a group of observations
    incrementally.
    """
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
        if self.n < 0:
            raise ValueError('Attempt to pull from an empty RunningVar')
        elif self.n == 0:
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
        if self.n == 0:
            #warnings.warn("merging two empty RVs.")
            self.m = 0
            self.s = 0
        else:
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
        if self.n < 0:
            raise ValueError('Attempt to split more obs from a RunningVar than it contains')
        elif self.n == 0:
            #warnings.warn('Splitting last observations from a RunningVar')
            self.m = 0
            self.s = 0
        else:
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


def sampleMultinomial(probs, n=1):
    """
    Sample from a multinomial distribution over given probabilities.
    Automatically normalizes probs to sum to 1, and returns a list of
    indices
    """
    # multinomial returns a trials-by-object matrix which has one
    # nonzero entry per row.  the second element returned by
    # nonzero() is the column (object) indices.
    return np.nonzero(np.random.multinomial(1, [float(x)/sum(probs) for x in probs], n))[1]


################################ TESTING STUFF ##########################################
def makeGaussianPhon(mean=0.0, var=1.0, n=20):
    """
    Sample some points from a gaussian distribution and return a Phon
    made from them as well as the sampled observations themselves

    Keyword args:
    mean=0.0
    var=1.0
    n=20
    """
    samp = np.random.normal(mean, var, n)
    p = Phon()
    for x in samp: p.add(x)
    return p, samp

def makeSomeWords():
    obs1 = np.random.normal(2, 1.0, 50)
    obs2 = np.random.normal(-2, 1.0, 50)
    words1 = zip(obs1[0:25], obs2[0:25])
    words2 = zip(obs2[25:50], obs1[25:50])
    words3 = [(x,) for x in np.random.normal(2, 1.0, 25)]
    words4 = [(x,) for x in np.random.normal(-2, 1.0, 25)]
    return (words1 + words2 + words3 + words4)

def makeWords1dToy(minPairs=False):
    A = list(np.random.normal(-5, 1.0, 400))
    B = list(np.random.normal(-1, 1.0, 200))
    C = list(np.random.normal(1, 1.0, 200))
    D = list(np.random.normal(5, 1.0, 400))
    if minPairs:
        wAB = [(A.pop(), B.pop()) for n in range(100)]
        wAC = [(A.pop(), C.pop()) for n in range(100)]
        wDB = [(D.pop(), B.pop()) for n in range(100)]
        wDC = [(D.pop(), C.pop()) for n in range(100)]
        words = wAB + wAC + wDB + wDC
    else:
        wAB = [(A.pop(), B.pop()) for n in range(200)]
        wDC = [(D.pop(), C.pop()) for n in range(200)]
        words = wAB + wDC
    wD = [(D.pop(), ) for n in range(100)]
    wADA = [(A.pop(), D.pop(), A.pop()) for n in range(100)]
    words += wD + wADA
    return words
