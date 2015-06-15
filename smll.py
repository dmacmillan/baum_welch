from baum_welch import *

S = ['3UTR','NULL']
N = len(S)
pi = Matrix.makeMatrix([[0.9,0.1]])

# Everything below is sequence specific
seq = 'AAAAAAAGT'
counts,O = countSeq(seq,6)

M = len(counts)
V = counts.keys()
A = Matrix(N,N)
A.mtx = [[0.9,0.1],[0.5,0.5]]
B = Matrix.makeRandom(N,M)
A.normalizeRows()
B.normalizeRows()
for i,v in enumerate(V):
    B.labels[v] = i
T = len(O)

for i in xrange(10):
    alpha,cs = calcAlpha(A,B,pi,O)
    beta = calcBeta(A,B,pi,O,cs)
    gamma = calcGamma(alpha,beta)
    xi = calcXi(alpha,beta,A,B,O)
    pi = calcNewPi(gamma,scale=False)
    A = calcNewA(xi,gamma)
    temp = calcNewB(gamma,O,V)
    temp.labels = B.labels
    B = temp
