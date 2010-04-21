import re, sys
from numpy import *
import lh2md


def readLemmas(fn='phondicts/vowels_only.dat'):
    infile = open(fn, 'r')
    # throw out header
    infile.readline()
    outdict = {}
    for line in infile:
        word, phons, vowels, freq = line.split(',')
        try: outdict[word] += [(phons, float(freq))]
        except KeyError: outdict[word] = [(phons, float(freq))]
    infile.close()
    return outdict

def makeVowelFreqDict(lemmadict, freqOnly=False):
    # consider anything with a stress marking a vowel
    vpat = re.compile(r'([A-Z]+)[012]')
    vowdict = {}
    for lemma in lemmadict:
        ndups = len(lemmadict[lemma])
        for phons, freq in lemmadict[lemma]:
            # get vowels, stripping off stress markers (may include later...)
            vowels = tuple(s.lower() for s in vpat.findall(phons))
            # adjust frequency for duplicate listings
            adjfreq = freq/ndups
            try: vowdict[vowels] += [(lemma, phons, adjfreq)]
            except KeyError: vowdict[vowels] = [(lemma, phons, adjfreq)]
    if freqOnly:
        for k, v in vowdict.iteritems():
            vowdict[k] = sum([x[2] for x in v])
    return vowdict

def freq(vowels, vdict):
    try: return sum([x[2] for x in vdict[vowels]])
    except KeyError: return None

def IPhODlexdict():
    """Generate a lexical dictionary compatible with makeWordsRand"""
    vlexdict = makeVowelFreqDict(readLemmas(), freqOnly=True)
    # we need to filter out words that have diphthongs
    goodvowels = cmu2hillDict.keys()
    gooddict = {}
    for lex in vlexdict:
        if all([p in goodvowels for p in lex]): gooddict[lex] = vlexdict[lex]
    return gooddict


def readHillenbrandVowels(fn='hillenbrand-data/vowdata.dat'):
    datlinepat = re.compile(r'^[mwbg][0-9]{2}')
    hillfile = open(fn, 'r')
    data = []
    for line in hillfile:
        # skip this line if it doesn't look like a data line
        if not datlinepat.search(line): continue
        stuff = line.split()
        name = stuff[0]       # which recording
        gender = name[0]
        vowel = name[3:5]
        dur = stuff[1]        # duration in ms
        ss = stuff[2:7]       # steady state f0, F1-F4
        s20 = stuff[7:10]     # F1-F3 at 20% duration
        s50 = stuff[10:13]    # F1-F3 at 50% duration
        s80 = stuff[13:16]    # F1-F3 at 80% duration
        # for the time being, just track the first two steady state formants (and ID info)
        data.append( (vowel, gender, array(ss[1:4], dtype=float)) )
    hillfile.close()
    return data

def crunchHill(hdata):
    vdict = {}
    for vowel, gender, F123 in hdata:
        vowel = hill2cmu(vowel)
        if gender=='b' or gender=='g': gender = 'c'
        try:
            vd = vdict[vowel]
            try: vd[gender].push(F123)
            except KeyError: vd[gender] = lh2md.MDRunningVar().push(F123)
        except KeyError:
            vdict[vowel] = {gender: lh2md.MDRunningVar().push(F123)}
    for vowel in vdict:
        All = lh2md.MDRunningVar()
        [All.merge(x) for x in vdict[vowel].values()]
        vdict[vowel]['all'] = All
    return vdict

def HillenbrandPhondict(gender='m', ndim=2):
    """Make the id->mvg() mapping to make words from..."""
    vdict = crunchHill(readHillenbrandVowels())
    phondict = {}
    if gender not in ['m', 'f', 'c', 'all']:
        raise ValueError("Gender must be 'm', 'f', 'c', or 'all'")
    for vowel, data in vdict.iteritems():
        mean = data[gender].mean[0:ndim]
        cov = data[gender].var()[0:ndim, 0:ndim]
        phondict[vowel] = lh2md.mvg(mean=mean, cov=cov)
    return phondict
        

hill2cmuDict = {'ae': 'ae', 'ah': 'aa', 
                'eh': 'eh', 'ei': 'ey',
                'ih': 'ih', 'iy': 'iy',
                'oa': 'ow', 'aw': 'ao',
                'uw': 'uw', 'oo': 'uh',
                'uh': 'ah', 'er': 'er'}
def hill2cmu(vowels):
    if isinstance(vowels, str): return hill2cmuDict[vowels]
    else: return [hill2cmuDict[x] for x in vowels]

cmu2hillDict = {'ae': 'ae', 'aa': 'ah', 
                'eh': 'eh', 'ey': 'ei',
                'ih': 'ih', 'iy': 'iy',
                'ow': 'oa', 'ao': 'aw',
                'uw': 'uw', 'uh': 'oo',
                'ah': 'uh', 'er': 'er'}
def cmu2hill(vowels):
    if isinstance(vowels, str): return cmu2hillDict[vowels]
    else: return [cmu2hillDict[x] for x in vowels]



################################################################################
## FINALLY make some words...
################################################################################

class DPRandomLexicon():
    """
    Represents a Dirichlet Process random lexicon, which provides observations
    through Chinese Restaurant Process sampling.
    """
    def __init__(self, base, concentration=10., phondict=HillenbrandPhondict()):
        self.base = base
        self.phons = phondict
        self.lexs = []
        self.lexcount = []
        self.con = float(concentration)

    def draw1(self, withLexs=False):
        counts = array(self.lexcount + [self.con])
        i = lh2md.sampleMultinomial(counts)
        if i == len(self.lexs):
            lex = self.base.draw()
            self.lexs.append(lex)
            self.lexcount.append(1)
        else:
            lex = self.lexs[i]
            self.lexcount[i] += 1
        samp = tuple(self.phons[ph].draw() for ph in lex)
        if withLexs: return (samp, lex)
        else: return samp

    def draw(self, n=None, withLexs=False):
        if n == None: return self.draw1(withLexs)
        else: return [self.draw1(withLexs) for i in range(n)]


class IPhODRandomLexicon(DPRandomLexicon):
    """
    Represents a random DP lexicon with a multinomial base distribution taken from
    the Irvine Phonetic Dictionary frequency counts.
    """
    def __init__(self, concentration=10., phondict=HillenbrandPhondict()):
        DPRandomLexicon.__init__(self,
                                 base=IPhODBaseDistribution(),
                                 concentration=concentration,
                                 phondict=phondict)


class IPhODBaseDistribution():
    """
    Represents a multinomial base distribution of lexemes with probabilities
    generated from frequencies in the Irvine Phonetic Dictionary
    """
    def __init__(self):
        lexdict = IPhODlexdict()
        self.lexs = lexdict.keys()
        totalprob = sum(lexdict.values())
        self.probs = array([x/totalprob for x in lexdict.values()])

    def draw(self, n=None):
        if n==None: return self.lexs[lh2md.sampleMultinomial(self.probs, 1)]
        else: return [self.lexs[i] for i in lh2md.sampleMultinomial(self.probs, n)]


class FeldmanRandomLexicon(DPRandomLexicon):
    """
    Represent a random lexicon based on the generative model from Feldman (2009)
    """
    def __init__(self, concentration=10., phondict=HillenbrandPhondict(), length=0.5):
        DPRandomLexicon.__init__(self,
                                 base=FeldmanBaseDistribution(phondict, length),
                                 concentration=concentration,
                                 phondict=phondict)


class FeldmanBaseDistribution():
    """
    Represents a base distribution defined by the generative model from Feldman
    (2009), where a length is randomly selected according to a geometric distribution
    and phones are independently assigned to each slot.
    """
    def __init__(self, phondict=HillenbrandPhondict(), length=0.5):
        self.phons = phondict.keys()
        self.phoncounts = [1.]*len(self.phons)
        self.len = length

    def drawPhon(self, n=None):
        if n == None: return self.phons[lh2md.sampleMultinomial(self.phoncounts)]
        else: return [self.phons[i] for i in lh2md.sampleMultinomial(self.phoncounts, n)]

    def draw1(self):
        length = random.geometric(self.len)
        return tuple(self.drawPhon(length))

    def draw(self, n=None):
        if n == None: return self.draw1()
        else: return [self.draw1() for n in range(n)]


hillHDPparams = {'nu': 1.001,
                 'mu': array([500., 1500.]),
                 's': array([[1., 0.], [0., 1.]]),
                 'alpha': 1,
                 'beta': 1,
                 'r': 5}
                 
