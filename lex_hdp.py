from numpy import *
from collections import deque

class Lex:

class Phon:

class Lexicon:
    def __init__(self, obs):
        self.lexs = deque( [] )
        self.obs = obs
    def iterate(self, ntimes):
        for n in range(ntimes):
            pass




class Phonicon:

class Observations:
    def __init__(self, obs):
        self.obs = obs
    def getObs(self, rows, cols):
        assert min(map(lambda r: obs[r]
        res = []
        for r in rows:
            for c in cols:
                res.append(self.obs[r][c]) 
