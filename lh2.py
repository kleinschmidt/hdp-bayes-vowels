class RunningVar:
    def __init__(self):
        self.n = 0
        self.m = 0.0
        self.s = 0.0
    
    def push(self, x):
        self.n += 1
        m_new = self.m + (x-self.m)/self.n
        self.s = ((self.n-1)*self.s + (x-self.m)*(x-m_new)) / self.n
        self.m = m_new
    
    def pull(self, x):
        self.n -= 1
        if (self.n==0):
            self.m = 0
            self.s = 0
        else:
            m_new = (self.m*(self.n+1) - x) / self.n
            self.s = (self.s*(self.n+1) - (x-self.m)*(x-m_new)) / self.n
            self.m = m_new
    
    def merge(self, other):
        n = self.n
        m = self.m
        s = self.s
    
        self.n = n + other.n
        self.m = (n*m + other.n*other.m) / self.n
        self.s = (n*s + other.n*other.s) / self.n + n*other.n/(n+other.n)**2 * (m-other.m)**2
    
    def split(self, other):
        n = self.n
        m = self.m
        s = self.s
    
        self.n -= other.n
        self.m = (n*m - other.n*other.m) / self.n
        self.s = s + other.n/self.n * (s-other.s) - other.n/n * (self.m-other.m)**2

