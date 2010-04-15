from scipy.special import multigammaln as lgammad
import numpy as np
import numbers, warnings, sys

class PhonDP:
    def __init__(self, params=None):
        if not params:
            self.params = { 'nu': 1.001, 'mu': 0.0, 's': 1.0, 'alpha': 1.0 }
        else:
            self.params = params
        self.children = None
        self.phons = [Phon(parent=self)]
        #self.priors = [Phon(parent=self, count=self.params['alpha']/self.params['s']) for n in range(self.params['s'])]
        # in a situation where the prior wasn't simple the above expression should
        # be used (so the prior probability is approximated by sampling
        # s items from the base distribution) but since the prior in this case is
        # analytically tractable (with an empty Phon) there's no need to bother
        self.priors = [Phon(parent=self, count=self.params['alpha'])]
    
    def newphon(self):
        p = Phon(self)
        self.phons.append(p)
        return p
    
    def iterate(self, n=1):
        for i in range(n):
            # iterate n times
            for seg in self.children.segments():
                # seg should be a Segment instance, check here
                if not isinstance(seg, Segment): raise TypeError('seg should be a Segment')
                # hold this segment out
                try:
                    lab = seg.holdout()
                except EmptyTable:
                    #...pruning it's label if necessary
                    print 'pruning Phon'
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
                    self.phons.append(seg.phon)
        return

    def prune(self, phon):
        """
        Remove a given phon from this DP (method exists pretty much only for
        symmetry's sake with LexDP
        """
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
            # append it to the list of lexs
            self.lexs.append(newlex)
            # add each word of the right length to it (and it's Phons...) and
            # append a properly labeled Word to self.words
            for w in wordsByLen[length]:
                newlex.add(w)
                self.words.append(Word(lex=newlex, obs=w))
            # initialize priors, too
            self.priors.extend(self.sampleBase(length=length, 
                                               n=self.params['r'], 
                                               count=self.params['beta']/self.params['r']))
    
    def iterate(self, n=1):
        for i in range(n):
            for word in self.words:
                # hold out this word
                try:
                    oldlex = word.holdout()
                except EmptyTable:
                    print 'pruning lex'
                    self.prune(oldlex)
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
                    # update lex counter and assign label to new Lex
                    newlex.id = self.lexcounter
                    self.lexcounter += 1
                    # record the new Lex
                    word.relabel(newlex)
                    self.lexs.append(newlex)
                    self.refreshPriors(length=len(word))
            # sweep through the parent PhonDP
            self.parent.iterate()
    
    def refreshPriors(self, length=None):
        """
        Resample Lexs of the specified length from the base distribution defined
        the parent PhonDP.
        """
        # get the indices to replace--either the indices w/ len==length or all
        # the priors if no length is specified
        indices = (i for i in range(len(self.priors))
                   if len(self.priors[i])==length or not length)
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
        for seg in lex:
            # throw in a sanity check, just for the hello of it.  if lex.count==0
            # then seg.count should == 0 for all segments.
            if seg.obs.n != lex.count:
                raise Exception("Segment (%d) and Lex (%d) counts disagree..." % 
                                (seg.obs.n, lex.count))
            # remove from parents, too
            try:
                seg.phon.remove(seg)
            except EmptyTable:
                self.parent.prune(seg.phon)
        self.lexs.remove(lex)
    
    def addLex(self, length):
        """Sample a new Lex and properly update it's parent Phons"""
        newlex = self.sampleBase(length, n=1, count=0).pop()
        for seg in newlex: seg.phon.add(seg)
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
    def __init__(self, parent=None, count=0, name="generic Phon"):
        self.obs = MDRunningVar()
        self.count = count
        self.name = name
        if parent==None:
            self.params = {'nu': 1.001, 'mu': 0.0, 's': 1.0}
        else:
            self.params = parent.params

    def __str__(self):
        return self.name

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
        if self.count == 0: raise EmptyTable('Prune this empty Phon!')
        if self.count < 0: raise ValueError('Phon count is less than 0...')

    def push(self, seg):
        """
        Push additional observations into this Phon, updating summary stats but
        not .count (use .add() instead)
        """
        if isinstance(seg, MDRunningVar):
            self.obs.merge(seg)
        elif isinstance(seg, np.ndarray):
            self.obs.push(seg)
        else:
            raise TypeError('seg needs to be ndarray or MDRunningVar to push to Phon (not %s)' % seg.__class__)

    def pull(self, seg):
        """
        Pull some observations from this Phon, updating summary stats but not
        .count (use .remove() instead)
        """
        if isinstance(seg, MDRunningVar):
            self.obs.split(seg)
        elif isinstance(seg, np.ndarray):
            self.obs.pull(seg)
        else:
            raise TypeError('seg needs to be ndarray or MDRunningVar to pull from Phon (not %s)' % seg.__class__)
    
    def lhood(self, seg, LOG=True):
        """
        Calcunlate the likelihood (defaults to log) of a given segment, which should
        be a MDRunningVar instance.
        """
        if not isinstance(seg, MDRunningVar):
            raise TypeError('seg must be a MDRunningVar instance (not %s)' % seg.__class__)
        n = seg.n
        nu = self._nu_n()
        mu = self._mu_n()
        s = self._s_n()

        d = float(max(self.obs.mean.shape))

        L = lgammad((nu+n)*0.5, d) + nu*0.5 * np.log(np.linalg.det(s)) - \
            lgammad(nu/2, d) - d*n*0.5 * np.log(np.pi) - d*0.5 * np.log((n+nu)/nu) - \
            (nu+n)*0.5 * np.log(np.linalg.det(s + seg.M2 + np.outer(seg.mean-mu, seg.mean-mu) * n*nu/(n+nu)))
        if LOG: return L
        else: return np.exp(L)

    def lhood1(self, seg, LOG=True):
        """
        Calculate the likelihood (defaults to log-lhood) for a single observation,
        which is either numeric or an ndarray.

        For likelihood of a group of observations, use .lhood().
        """
        return self.lhood(MDRunningVar(n=1, mean=seg, M2=0.), LOG)
    
    def _nu_n(self):
        return self.params['nu'] + self.obs.n
    
    def _mu_n(self):
        return (self.params['mu']*self.params['nu'] + self.obs.mean*self.obs.n) / (self.params['nu']+self.obs.n)
    
    def _s_n(self):
        return (self.params['s'] + self.obs.M2 + 
                self.params['nu']*self.obs.n/(self.params['nu']+self.obs.n) * \
                    np.outer(self.obs.mean-self.params['mu'], self.obs.mean-self.params['mu']))
    
    def test(self):
        # this checks out with the R implementation in lh2.r
        x = [Segment(MDRunningVar(n=1,mean=float(x),M2=0.), self) for x in range(20)]
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
        ret = 'Lex: ' + '\'untitled\'' if not self.id else str(self.id)
        for seg in self:
            ret += '\n' + seg.__str__()
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
        if self.count == 0: raise EmptyTable('Prune this empty Lex!')
        if self.count < 0: raise ValueError('Lex count should not be less than 0...')
    
    def holdout(self, word):
        self.remove(word)
    
    def lhood(self, word):
        if len(self) != len(word):
            return 0
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
        self.lex.remove(self.obs)
        l = self.lex
        self.lex = None
        return l


class Segment:
    """
    Simple container class, to bind together a single observation and a label.

    Fields:
    obs (should be a RunningVar instance)
    phon (should be a Phon instance)
    """
    def __init__(self, obs, phon):
        if not isinstance(obs, MDRunningVar):
            raise TypeError('Segment obs should be a MDRunningVar instance')
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


class MDRunningVar():
    def __init__(self, n=0, mean=0., M2=0.):
        # count of observations
        self.n = n
        # mean of observations: make 1-D array of proper size
        self.mean = np.array(mean).reshape(-1)
        # second cumulant (\sum_i (x_i-\bar{x})^2)
        if M2==0.:
            self.M2 = np.zeros(self.mean.shape*2)
        else:
            try:
                self.M2 = np.array(M2).reshape(self.mean.shape * 2)
            except ValueError:
                # this should catch cases where the dimensions of the mean and covariance
                # matrices mismatch
                raise ValueError('Mean and covariance matrices must have congruent sizes')
            
    def __str__(self):
        return "mrv(n=%s, m=%s, s=%s)" % (self.n, self.mean, self.var())

    def __repr__(self):
        return "mrv(n=%s, m=%s, s=%s)" % (repr(self.n), repr(self.mean), repr(self.var()))
        
    def var(self):
        """Convert second cumulant (which is tracked incrementally) into variance"""
        if self.n > 0: return self.M2 / self.n
        else: return np.zeros(self.M2.shape)

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


def sampleMultinomial(probs, n=1):
    """
    Sample from a multinomial distribution over given probabilities.
    Automatically normalizes probs to sum to 1, and returns a list of
    indices
    """
    # multinomial returns a trials-by-object matrix which has one
    # nonzero entry per row.  the second element returned by
    # nonzero() is the column (object) indices.
    return np.nonzero(np.random.multinomial(1, [x/sum(probs) for x in probs], n))[1]


def debug_methodInfo(func):
    def wrapper(self, *args, **kwargs):
        sys.stdout.write( '%s.%s ' % (self.__class__.__name__, func.__name__) )
        func(self, *args, **kwargs)
        sys.stdout.write( ' (done)\n' )
    return wrapper

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
    return (words1 + words2)
