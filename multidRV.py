import numpy as np

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

class MDRunningVar():
    def __init__(self):
        self.n = 0
        self.m = 0.0
        self.s = 0.0

    def push(self, x):
        self.n += 1
        delta = x-self.m
        self.m += delta/self.n
        self.s += (self.s - np.outer(delta, (x-self.m))) / self.n

    def pull(self, x):
        
