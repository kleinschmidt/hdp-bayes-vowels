import lh2md, lexgen
import numpy as np

####################INFO THEORY HELPER FUNCTIONS################################

def mylog2(x):
    """Calculate base-2 log, defining log(0) = 0"""
    X = array(x)
    ret = zeros(X.shape)
    ret[X!=0] = log(X[X!=0])/log(2)
    return ret

def mutinf(JPD):
    """Calculate the mutual information of a joint probability distribution"""
    jpd = array(JPD, dtype=float)
    jpd /= jpd.sum()
    return sum(jpd* (mylog2(jpd) - mylog2(outer(jpd.sum(axis=1), jpd.sum(axis=0)))))

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
    print '\n'
    print 'real  learned'
    print '     ' + ' '.join(['%4s' % lab for lab in ldpphons])
    for real in realphons:
        print ('%4s ' % real) + ' '.join(['%4d' % paircounts[(real,lab)] for lab in ldpphons])
    # put the output into a big array: confarray[i,j] = how often realphon[i] was labeled ldpphons[j]
    confarray = np.array([paircounts[a,b] for a in realphons for b in ldpphons]).reshape(len(realphons),len(ldpphons))
    return confarray

def ldpMutualInfo(ldp, lexlabels):
    """Calculate the mutual information of the actual and learned phon labels"""
    return mutinfo(confusionTable(ldp,lexlabels))

###########################SIMULATION FUNCTIONS#################################
