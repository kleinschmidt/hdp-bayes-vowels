import lh2md, lexgen
import numpy as np

####################INFO THEORY HELPER FUNCTIONS################################

def mylog2(x):
    """Calculate base-2 log, defining log(0) = 0"""
    X = np.array(x)
    ret = np.zeros(X.shape)
    ret[X!=0] = np.log(X[X!=0])/np.log(2)
    return ret

def mutinf(JPD):
    """Calculate the mutual information of a joint probability distribution"""
    jpd = np.array(JPD, dtype=float)
    jpd /= jpd.sum()
    return (jpd* (mylog2(jpd) - mylog2(np.outer(jpd.sum(axis=1), jpd.sum(axis=0))))).sum()

#################LEXDP ANALYSIS FUNCTIONS#######################################

def confusionTable(ldp, lexlabels):
    """
    Generate a phon-by-phon "confusion" table from an LexDP and a list of real
    labels
    """
    ldpphons = [p.id for p in ldp.parent.phons]
    realphons = set(p for l in lexlabels for p in l)
    keys = [(a,b) for a in realphons for b in ldpphons]
    paircounts = dict((k, 0) for k in keys)
    realflat = [a for l in lexlabels for a in l]
    ldpflat = [seg.phon.id for word in ldp.words for seg in word.lex]
    for pair in zip(realflat, ldpflat):
        paircounts[pair] += 1
    # print output:
    print 'real  learned'
    print '     ' + ' '.join(['%4s' % lab for lab in ldpphons])
    for real in realphons:
        print ('%4s ' % real) + ' '.join(['%4d' % paircounts[(real,lab)] for lab in ldpphons])
    # put the output into a big array: confarray[i,j] = how often realphon[i] was labeled ldpphons[j]
    confarray = np.array([paircounts[a,b] for a in realphons for b in ldpphons]).reshape(len(realphons),len(ldpphons))
    return confarray

def phonMutualInfo(ldp, lexlabels):
    """Calculate the mutual information of the actual and learned phon labels"""
    return mutinf(confusionTable(ldp,lexlabels))

###########################SIMULATION FUNCTIONS#################################

def simulate(niter=100, printEvery=-1, numWords=5000, randLex=lexgen.FeldmanRandomLexicon):
    """
    Run the Feldman lexical-phonetic learner for specified number of iterations
    with running analysis of phon confusions, mutual info, etc.  Keyword arguemnt
    randLex specifies which class to use to generate observations.  It needs to
    provide a method .draw(n, withLexs) method that returns tuples of observations
    and lexical labels (or only observations if withLexs==False) and defaults to
    lexgen.CRPRandomLexicon which implements Feldman's generative model.

    TODO: implement precision and recall calculation...
    """
    # instantiate randLex and draw observations and labels
    print 'Simulating LexDP for %d iterations on %d words from %s' % \
        (niter, numWords, randLex.__name__)
    print '...drawing samples'
    lexicon = randLex()
    samp, lexs = zip(*lexicon.draw(numWords, withLexs=True))
    # create and initialize a LexDP
    ldp = lh2md.LexDP(samp, params=lexgen.hillHDPparams)
    for n in range(niter):
        print 'Iteration', n+1
        ldp.iterate()
        if mod(n+1,printEvery)==1: print ldp
        print 'Phon confusion table:'
        print 'Phon mut. info.: %.4f' % phonMutualInfo(ldp, lexs)
    return ldp
