from baum_welch import *

S = ['3UTR','NULL']
N = len(S)
pi = Matrix.makeMatrix([[0.9,0.1]])

# Everything below is sequence specific
seq = 'AAAAAAAGT'
counts,O = countSeq(seq,6)

M = len(counts)
V = counts.keys()
A = Matrix.makeRandom(N,N)
B = Matrix.makeRandom(N,M)
A.normalizeRows()
B.normalizeRows()
for i,v in enumerate(V):
    B.labels[v] = i
T = len(O)

print 'A:'
print A
print
print 'B:'
print B
print

alpha = calcAlpha(A,B,pi,O,scale=False)[0]
print 'alpha (not scaled):'
print alpha
print

salpha,cs = calcAlpha(A,B,pi,O)
print 'alpha (scaled):'
print salpha
print

beta = calcBeta(A,B,pi,O)
print 'beta (not scaled):'
print beta
print

sbeta = calcBeta(A,B,pi,O,cs)
print 'beta (scaled):'
print sbeta
print

# This will be the same whether we use
# the scaled alpha/beta or not
gamma = calcGamma(salpha,sbeta)
print 'gamma:'
print gamma
print

xi = calcXi(salpha,sbeta,A,B,O)
print 'xi:'
print xi
print

newpi = calcNewPi(gamma,scale=False)
print 'newpi:'
print newpi
print

newa = calcNewA(xi,gamma)
print 'newa:'
print newa
print

newb = calcNewB(gamma,B,O,counts)
print 'newb:'
print newb
print
