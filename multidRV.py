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
        self.n = 0          # count of observations
        self.mean = 0.0     # mean of observations
        self.M2 = 0.0       # second cumulant (\sum_i (x_i-\bar{x})^2)

    def __str__(self):
        return "mrv(n=%s, m=%s, s=%s)" % (self.n, self.m, self.var())

    def __repr__(self):
        return "mrv(n=%s, m=%s, s=%s)" % (repr(self.n), repr(self.m), repr(self.var()))
        
    def var(self):
        """Convert second cumulant (which is tracked incrementally) into variance"""
        if self.n > 0: return self.M2 / self.n
        else: return 0

    def push(self, x):
        delta = x - self.mean
        self.n += 1
        self.mean += delta/self.n
        self.M2 += np.outer(delta, (x-self.mean))

    def pull(self, x):
        delta = x - self.mean
        self.n -= 1
        self.mean -= delta/self.n
        self.M2 -= np.outer((x-self.mean), delta)

    def merge(self, other):
        delta = other.mean - self.mean
        M2delta = np.outer(delta, delta) * self.n*other.n / (self.n+other.n)
        self.n += other.n
        self.mean += delta*other.n/self.n
        self.M2 += other.M2 + M2delta
