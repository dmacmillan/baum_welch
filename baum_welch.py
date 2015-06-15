import random,time,sys

class HMM:

    def __init__(self,A,B,pi):
        self.A = A
        self.B = B
        self.pi = pi

class MatrixError(Exception):
    pass

class Matrix:

    def __init__(self,m,n,z=None,init=True,labels={}):
        self.mtx = []
        if init:
            if not z:
                self.mtx = [[0]*n for x in xrange(m)]
            else:
                self.mtx = [[[0 for i in xrange(z)] for j in xrange(n)] for k in xrange(m)]
        self.m = m
        self.n = n
        self.z = z
        self.labels = labels

    def __str__(self):
        s='\n'.join([' '.join([str(item) for item in row]) for row in self.mtx])
        return s

    def __mul__(self,mtx):
        if (self.n != mtx.m):
            raise MatrixError, 'Matrices cannot be multiplied!'
        pass

    def col(self,c):
        return [x[c] for x in self.mtx]

    def setCol(self,c,v):
        if len(v) != self.m:
            raise MatrixError, 'Vector must be same size as n for an mxn matrix'
        for row in xrange(self.m):
            self.mtx[row][c] = v[row]

    def row(self,r):
        return self.mtx[r]

    def xy(self,x,y):
        return self.mtx[x][y]
    # This will normalize the rows such that all elements
    # in each row will add to 1
    def normalizeRows(self):
        for i in xrange(len(self.mtx)):
            multiplier = float(1) / sum(self.mtx[i])
            self.mtx[i] = [x*multiplier for x in self.mtx[i]]
    
    @classmethod
    def dot(cls,v1,v2):
        return sum(p*q for p,q in zip(v1,v2))

    #@classmethod
    #def multiply(cls,*args):
        
    # Make an m row by n column (mxn) matrix
    # of random values between low and high 
    @classmethod
    def makeRandom(cls,m,n,labels={}):
        obj = Matrix(m,n,init=False,labels=labels)
        for x in xrange(m):
            obj.mtx.append([random.random() for i in xrange(obj.n)])
        return obj

    @classmethod
    def makeMatrix(cls,mtx,labels={}):
        m = len(mtx)
        n = len(mtx[0])
        if any([len(row) != n for row in  mtx[1:]]):
            raise MatrixError, "Inconsistent row length"
        mat = Matrix(m,n,init=False)
        mat.mtx = mtx
        mat.labels = labels
        return mat

# For an alphabet A (string), return a list
# of all possible strings of length x (int)
def generateAllXmers(A,x):
    c = []
    for i in xrange(x-1):
        c = [x+y for x in A for y in c or A]
    return c

#sys.stdout.write("Normalizing rows... ")
#s = time.time()
#sys.stdout.write("Done. Time: {}\n".format(time.time()-s))

def dot3(x,y,z):
    return sum([a*b*c for a,b,c in zip(x,y,z)])

def calcAlpha(A,B,pi,obs,scale=True):
    alpha = Matrix(B.m,len(obs))
    cs = []
    #print 'Created alpha = {}x{} matrix'.format(alpha.m,alpha.n)
    for t in xrange(alpha.n):
        #print '  Time: {}'.format(t)
        for s in xrange(alpha.m):
            #print '  State: {}'.format(s)
            if (t == 0):
                alpha.mtx[s][0] = pi.xy(0,s)*B.xy(s,B.labels[obs[0]])
                #print '    alpha[{}][0] = {}'.format(s,alpha.xy(s,t))
            else:
                alpha.mtx[s][t] = B.xy(s,B.labels[obs[t]])*Matrix.dot(alpha.col(t-1),A.col(s))
                #print '    alpha[{}][{}] = {}'.format(s,t,alpha.xy(s,t))
        #for s in xrange(alpha.m):
        #    c = 1 / sum(alpha.col(t))
        #    alpha.mtx[s][t] *= c
        if scale and (t > 0):
            c = 1 / sum(alpha.col(t))
            cs.append(c)
            for s2 in xrange(alpha.m):
                #print '  Scaling:'
                #print '    alpha[{}][{}] = {} / sum({})'.format(s,t,alpha.xy(s,t),alpha.col(t))
                alpha.mtx[s2][t] *= c
                #print '    alpha[{}][{}] = {}'.format(s,t,alpha.xy(s,t))
        #raw_input('*')
    return [alpha,cs]

def calcBeta(A,B,pi,obs,cs=None):
    beta = Matrix(B.m,len(obs))
    #print 'Created beta = {}x{} matrix'.format(beta.m,beta.n)
    # Set for t = last element
    for s in xrange(beta.m):
        beta.mtx[s][beta.n-1] = 1
        #print '  Set beta[{}][{}] = 1'.format(s,beta.n-1)
    # Scale these
    #if cs:
    #    for s in xrange(beta.m):
    #        beta.mtx[s][beta.n-1] *= cs[beta.n-1]
            #print '  Scaled beta[{}][{}] *= {} = {}'.format(s,beta.n-1,cs[beta.n-1],beta.xy(s,beta.n-1))
    for t in xrange(beta.n-2,-1,-1):
        #print '  Time: {}'.format(t)
        for s in xrange(beta.m):
            #print '    State: {}'.format(s)
            beta.mtx[s][t] = dot3(beta.col(t+1),A.mtx[s][:],B.col(B.labels[obs[t+1]]))
            #print '      beta[{}][{}] = {}'.format(s,t,beta.xy(s,t))
        # Scale remaining
        if cs:
            for s2 in xrange(beta.m):
                beta.mtx[s2][t] *= cs[t]
        #raw_input('*')
    return beta

def calcGamma(alpha,beta):
    gamma = Matrix(alpha.m,alpha.n)
    for s in xrange(gamma.m):
        for t in xrange(gamma.n):
            num = alpha.xy(s,t)*beta.xy(s,t)
            denom = 0
            for s2 in xrange(gamma.m):
                denom += (alpha.xy(s2,t)*beta.xy(s2,t))
            gamma.mtx[s][t] = num/denom
    return gamma

def calcXi(alpha,beta,A,B,O):
    xi = Matrix(len(O)-1,B.m,B.m)
    for t in xrange(xi.m):
        for i in xrange(xi.n):
            for j in xrange(xi.z):
                #print 'Setting xi[{}][{}][{}]'.format(t,i,j)
                num = alpha.xy(i,t)*A.xy(i,j)*B.xy(j,B.labels[O[t+1]])*beta.xy(j,t+1)
                #print '  num = {}*{}*{}*{}'.format(alpha.xy(i,t),A.xy(i,j),B.xy(j,B.labels[O[t+1]]),beta.xy(j,t+1))
                denom = 0
                for p in xrange(xi.n):
                    for q in xrange(xi.z):
                        #print '  denom += (({})({})({})({}))'.format(alpha.xy(p,t),A.xy(p,q),B.xy(q,B.labels[O[t+1]]),beta.xy(q,t+1))
                        denom += (alpha.xy(p,t)*A.xy(p,q)*B.xy(q,B.labels[O[t+1]])*beta.xy(q,t+1))
                        #print '    denom += (alpha[{}][{}])(a[{}][{}])(b[{}][{}])(beta[{}][{}]) = {}'.format(p,t,p,q,q,B.labels[O[t+1]],q,t+1,denom)
                        #print '    denom = {}'.format(denom)
                xi.mtx[t][i][j] = num/denom
                #raw_input('*')
    return xi

#def calcXi(alpha,beta,A,B,obs):
#    xi = Matrix(B.m,B.m,len(obs)-1)
#    for s1 in xrange(xi.m):
#        for s2 in xrange(xi.n):
#            for t in xrange(xi.z):
#                num = alpha.xy(s1,t)*A.xy(s1,s2)*beta.xy(s2,t+1)*B.xy(s2,B.labels[obs[t]])
#                denom = 0
#                for k in xrange(xi.m):
#                    denom += alpha.xy(k,t)*beta.xy(k,t)
#                xi.mtx[s1][s2][t] = num/denom
#    return xi

def calcNewPi(gamma,scale=True):
    newpi = Matrix(1,gamma.m)
    for s in xrange(gamma.m):
        newpi.mtx[0][s] = gamma.xy(s,0)
    if scale:
        c = 1 / sum(newpi.mtx[0])
        for s in xrange(newpi.n):
            newpi.mtx[0][s] *= c
    return newpi

def calcNewA(xi,gamma):
    newa = Matrix(gamma.m,gamma.m)
    for i in xrange(newa.m):
        #print 'i: {}'.format(i)
        denom = sum(gamma.mtx[i][:-1])
        #print 'denom = sum({}) = {}'.format(gamma.mtx[i],denom)
        num = 0
        for j in xrange(newa.n):
            #print '  j: {}'.format(j)
            #for x in xi.mtx:
                #num += x[i][j]
            #print 'num = sum([x[{}] for x in xi.col({})])'.format(j,i)
            #print '=sum({})'.format([x[j] for x in xi.col(i)])
            num = sum([x[j] for x in xi.col(i)])
            newa.mtx[i][j] = num/denom
            #print 'newa[{}][{}] = {}/{}'.format(i,j,num,denom)
            #print 'newa[{}][{}] = {}/{} = {}'.format(i,j,num,denom,(num/denom))
    return newa

def calcNewB(gamma,O,V):
    vlen = len(V)
    newb = Matrix(gamma.m,vlen)
    for t in xrange(newb.n):
        for s in xrange(newb.m):
            newb.mtx[s][t] = sum([gamma.mtx[s][i] for i in xrange(gamma.n) if O[i] == V[t]]) / sum(gamma.mtx[s])
    return newb
    
# Compute the probability of
# an hmm
def Phmm(alpha,T):
    return sum(alpha.col(T-1))

def compAlpha(A,B,pi,O):
    alpha = {}
    for o in O.mtx:
        if o in alpha:
            continue
        else:
            alpha[o] = Matrix(A.m,len(o))
        for t in xrange(len(o)):
            for s in xrange(A.m):
                if (t == 0):
                    alpha[o].mtx[s][t] = pi.xy(s,0)*B.xy(s,O.labels[o[t]])
                else:
                    alpha[o].mtx[s][t] = B.xy(s,O.labels[o[t]])*Matrix.dot(alpha[o].col(t-1),A.col(s))
    return alpha

def compBeta(A,B,pi,O):
    beta = {}
    for o in O.mtx:
        olen = len(o)
        if o in beta:
            continue
        else:
            beta[o] = Matrix(A.m,olen)
        for s in xrange(A.m):
            beta[o].mtx[s][olen-1] = 1
        for t in xrange(olen-2,-1,-1):
            for s in xrange(A.m):
                beta[o].mtx[s][t] = dot3(beta[o].col(t+1),A.mtx[s][:],B.col(O.labels[o[t+1]]))
    return beta

def compGamma(alpha,beta,A,B,O):
    gamma = {}
    for o in O.mtx:
        denom = sum(alpha[o].col(alpha[o].n-1))
        olen = len(o)
        if o in gamma:
            continue
        else:
            gamma[o] = Matrix(olen-1,A.m,A.m)
        for t in xrange(gamma[o].m):
            for s1 in xrange(gamma[o].n):
                for s2 in xrange(gamma[o].z):
                    num = alpha[o].xy(s1,t)*A.xy(s1,s2)*beta[o].xy(s2,t+1)*B.xy(s2,O.labels[o[t+1]])
                    gamma[o].mtx[t][s1][s2] = num/denom
    return gamma

def compSigma(alpha,gamma,O):
    sigma = {}
    for o in O.mtx:
        olen = len(o)
        if o in sigma:
            continue
        else:
            sigma[o] = Matrix(alpha[o].m,alpha[o].n)
        for t in xrange(sigma[o].n-1):
            for s in xrange(sigma[o].m):
                sigma[o].mtx[s][t] = sum(gamma[o].mtx[t][s])
        for s in xrange(sigma[o].m):
            sigma[o].mtx[s][-1] = alpha[o].xy(s,alpha[o].n-1)/sum(alpha[o].col(alpha[o].n-1))
    return sigma

def countObs(O):
    d = {}
    for o in O.mtx:
        if o in d:
            d[o] += 1
        else:
            d[o] = 1
    return d

# return the counts of every sub-sequence
# of length k within sequence s as well
# as the actual subsequences
def countSeq(s,k):
    res = {}
    l = []
    for i in xrange(len(s)-k+1):
        l.append(s[i:i+k])
        if s[i:i+k] not in res:
            res[s[i:i+k]] = 1
        else:
            res[s[i:i+k]] += 1
    return [res,l]

def compI(sigma,counts):
    res = []
    for s in xrange(sigma.itervalues().next().m):
        i = 0
        for o in counts:
            i += (sigma[o].xy(s,0)*counts[o])
        res.append(i)
    return res

def compNewA(gamma,counts):
    numstates = gamma.itervalues().next().n
    newA = Matrix(numstates,numstates)
    for s0 in xrange(numstates):
        for s1 in xrange(numstates):
            for s2 in xrange(numstates):
                val = 0
                for o in counts:
                    i = 0
                    for t in xrange(gamma[o].m):
                        i += gamma[o].mtx[t][s1][s2]
                    i *= counts[o]
                    val += i
                newA.mtx[s1][s2] = val
    res = Matrix(newA.m,newA.n)
    for s1 in xrange(res.m):
        for s2 in xrange(res.n):
            if (s2 != res.n-1):
                res.mtx[s1][s2] = newA.xy(s1,s2)/sum(newA.mtx[s1])
            else:
                res.mtx[s1][s2] = 1 - sum(res.mtx[s1])
    return res

def compNewB(sigma,counts,B,O):
    newB = Matrix(B.m,B.n)
    for s in xrange(newB.m):
        for y in xrange(newB.n-1):
            d = {}
            for o in counts:
                if o not in d:
                    d[o] = {'in': 0, 'out': 0}
                for t in xrange(len(o)):
                    if (O.labels[o[t]] == y):
                        d[o]['in'] += sigma[o].xy(s,t)
                    else:
                        d[o]['out'] += sigma[o].xy(s,t)
                d[o]['in'] *= counts[o]
                d[o]['out'] *= counts[o]
            allin = sum([d[x]['in'] for x in d])
            allout = sum([d[x]['out'] for x in d])
            newB.mtx[s][y] =  allin / (allin + allout) 
        newB.mtx[s][-1] = 1 - sum(newB.mtx[s])
    return newB
