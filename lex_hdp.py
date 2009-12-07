from numpy import *

class Phones:
    def __init__(self):
        self.phones = []

    def addPhone(self, lex, pos):
        self.phones.append( (lex,pos) )

    # remove empty phoneme entries
    def clean(self):
        for p in self.phones:
            if len(p) == 0:
                self.phones.remove(p)

class Lexicon:
    def __init__(self):
        self.lex = []
        self.obs = array( [] )
