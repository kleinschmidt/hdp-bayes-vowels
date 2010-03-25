class RunningVar:
    def __init__(self, n=0, m=0.0, s=0.0):
        self.n = n
        self.m = m
        self.s = s

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
        if (self.n==0):
            self.m = 0
            self.s = 0
        else:
            m_new = (self.m*(self.n+1) - x) / self.n
            self.s = (self.s*(self.n+1) - (x-self.m)*(x-m_new)) / self.n
            self.m = m_new
    
    def merge(self, other):
        """merge the observations from another rv into this one"""
        nx = self.n
        mx = self.m
        sx = self.s

        ny = other.n
        my = other.m
        sy = other.s

        self.n = nx + ny
        self.m = (nx*mx + ny*my)/self.n
        self.s = (nx*sx + ny*sy)/self.n + (mx-my)**2 * nx*ny/(nx+ny)**2
    
    def split(self, other):
        """remove multiple observations from this rv"""
        n = self.n
        m = self.m
        s = self.s
        
        self.n -= other.n
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

