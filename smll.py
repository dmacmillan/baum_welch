from baum_welch import *

S = ['3UTR','NULL']
N = len(S)
#pi = Matrix.makeMatrix([[0.9,0.1]])

# Read in 3utrs
f = open('realsubset','r')
lines = f.readlines()
lines = [x.strip().split('\t') for x in lines]
lines = [x[2] for x in lines]
Os = lines[:]
lines = None

# Generate all possible 6-mers
V = generateAllXmers('ATCG',6)

# Number of possible observations
M = len(V)

finalhmm = HMM
finalhmm.A = Matrix.makeMatrix([[0.9,0.1],[0.6,0.4]])
finalhmm.B = Matrix.makeRandom(N,M)
for i in xrange(M):
        finalhmm.B.labels[V[i]] = i
finalhmm.pi = Matrix.makeMatrix([[0.9,0.1]])
#finalhmm.A.normalizeRows()
finalhmm.B.normalizeRows()

hmms = []
# Iterate over two sequences to start
#for k,seq in enumerate(Os[:2]):
for k,seq in enumerate(Os):
    counts,O = countSeq(seq,6)
    alpha,cs = calcAlpha(finalhmm.A,finalhmm.B,finalhmm.pi,O)
    beta = calcBeta(finalhmm.A,finalhmm.B,finalhmm.pi,O,cs)
    xi = calcXi(alpha,beta,finalhmm.A,finalhmm.B,O)
    gamma = calcGamma(alpha,beta)
    newhmm = HMM(None,None,None,init=False)
    newhmm.cs = cs
    newhmm.O = O
    newhmm.counts = counts
    newhmm.alpha = alpha
    newhmm.beta = beta
    newhmm.xi = xi
    newhmm.gamma = gamma
    hmms.append(newhmm)
    
for i in xrange(10):
        for hmm in hmms:
            hmm.alpha,hmm.cs = calcAlpha(finalhmm.A,finalhmm.B,finalhmm.pi,hmm.O)
            hmm.beta = calcBeta(finalhmm.A,finalhmm.B,finalhmm.pi,hmm.O,hmm.cs)
            hmm.xi = calcXi(hmm.alpha,hmm.beta,finalhmm.A,finalhmm.B,hmm.O)
            hmm.gamma = calcGamma(hmm.alpha,hmm.beta)
        finalhmm.A = multiCalcNewA(finalhmm,hmms)
        finalhmm.B = multiCalcNewB(finalhmm,hmms,V)
    
#T = len(O)

#for i in xrange(1):
#    alpha,cs = calcAlpha(A,B,pi,O)
#    beta = calcBeta(A,B,pi,O,cs)
#    gamma = calcGamma(alpha,beta)
#    xi = calcXi(alpha,beta,A,B,O)
#    pi = calcNewPi(gamma,scale=False)
#    A = calcNewA(xi,gamma)
#    temp = calcNewB(gamma,O,V)
#    temp.labels = B.labels
#    B = temp
