#!/usr/bin/env python

import sys
import numpy as np
import networkx as nx
import sklearn as sk
from sklearn.cluster import KMeans
import multiprocessing as mp
from scipy import spatial

class dist:
    def __init__(self, X):
        self.X = X
        (self.n, self.d) = X.shape
        sz = self.n*(self.n-1)/2
        Z = np.dot(X, X.T)
        self.P = np.zeros(sz)
        #for i in range(n-1):
        #    for j in range(i+1,n):
        #        r2sum = Z[i,i]+Z[j,j]-2*Z[i,j] 
        #        self.P[self.idx(i,j)] = np.sqrt(r2sum)
        ones = np.ones((n,1))
        diag = np.diag(Z)
        diag.shape = (n,1)
        L = np.dot(diag, ones.T)
        R = np.dot(ones, diag.T)
        self.P = np.sqrt(L + R - 2*Z) #proximity matrix
        self.D = np.diagflat(np.sum(self.P, axis=0)) #degree matrix
        self.L = self.D - self.P #laplacian matrix
    
    def idx(self, i,j):    
        idxx = ((self.n-1)*i + j - (i+1)*i/2 -1)
        return idxx
    
    def get_dist(self, i,j):
        if j>i:
            return self.P[self.idx(i,j)]
        else:
            return None

    def cluster_stats(self, C, k, nc, out=True):
        self.k = k
        self.C = C
        self.nc = nc #cluster sizes
        self.c = np.zeros((n,k)) #cluster indicator vector
        for i in range(n):
            self.c[i,C[i]] = 1 

        self.W = 0.0
        self.B = 0.0
        self.NW = 0
        self.NB = 0
        self.Wi = np.zeros(k)
        self.Bi = np.zeros(k)
        self.voli = np.zeros(k)
        for i in range(k):
            self.NW += nc[i]*(nc[i]-1)/2
            ci = self.c[:,i]
            ci.shape = (n,1)
            self.Wi[i] = np.dot(ci.T, np.dot(self.P, ci))
            self.Bi[i] = np.dot(ci.T, np.dot(self.L, ci))
            self.voli[i] = np.dot(ci.T, np.dot(self.D, ci))
            self.W += self.Wi[i]
            self.B += self.Bi[i]
            for j in range(i+1,k):
                self.NB += nc[i]*nc[j]
        
        #verify W,B
        WW = np.zeros((self.k, self.k))
        for i in range(self.n):
            for j in range(self.n):
                WW[C[i], C[j]] += self.P[i,j]

        if out:
            print "IN Cluster Stats"
            print WW
            print self.Wi, self.W
            print self.Bi, self.B
            print self.voli
        #print WW, self.W, 2*WW
        #print BB, self.B, 2*BB
        #assert(2*WW-self.W < 0.001)
        #assert(2*BB-self.B < 0.001)

        #cluster means
        self.mu = np.zeros((self.k, self.d))
        for i in range(n):
            cl = self.C[i]
            self.mu[cl,:] += X[i,:]
        self.mean = np.sum(self.mu, axis=0)
        for i in range(k):
            self.mu[i,:] /= self.nc[i]
        self.mean /= self.n
        if out:
            print "cluster means", self.mu
            print "total mean", self.mean

        #scatter matrices
        SW = np.zeros((self.d, self.d))
        SB = np.zeros((self.d, self.d))
        ST = np.zeros((self.d, self.d))
        
        for i in range(self.n):
            cl = self.C[i]
            dd = self.X[i,:]-self.mu[cl,:]
            dd.shape  = (self.d,1)
            SW += np.dot(dd, dd.T)
            dd = self.X[i,:]-self.mean
            dd.shape  = (self.d,1)
            ST += np.dot(dd, dd.T)
        for i in range(self.k):
            dd = self.mu[i,:] - self.mean
            dd.shape  = (self.d,1)
            SB += self.nc[i]*np.dot(dd, dd.T)

        self.SW = SW
        self.SB = SB
        self.ST = ST
        if out:
            print "SCATTERS"
            print SW
            print SB
            print ST

    #return the sum of the m min and max distances
    def sum_min_max(self, m):
        TU = []
        for i in range(self.n-1):
            for j in range(i+1, self.n):
                TU.append(self.P[i,j])
        
        sPary = np.sort(TU, axis=None)
        minM = 0.0
        maxM = 0.0
        for i in range(int(m)):
            minM += sPary[i] 
            maxM += sPary[-(i+1)]
        return (minM, maxM)
            
    #silhouette_score for i-th point
    def silhouette(self, i, out=False):
        ss = np.zeros(self.k)
        for j in range(self.n):
            ss[self.C[j]] += self.P[i,j]
        for j in range(self.k):
            if j == self.C[i]:
                ss[j] /= (self.nc[j]-1)
            else:
                ss[j] /= self.nc[j]

        wi = ss[self.C[i]]
        if self.k > 1:
            bi = np.min([ss[x] for x in range(self.k) if x != self.C[i]])
        else: 
            bi = 0.0
        if i == 0 and out: #just print out for one point
            print "Sil", i, self.C[i], ss, wi, bi, (bi-wi)/np.max([bi, wi])
        return ((bi-wi)/np.max([bi, wi]))

    def dunn(self):
        Bmin = sys.maxint
        Wmax = 0.0
        maxx = np.zeros((self.k, self.k))
        minn = np.ones((self.k, self.k)) * sys.maxint
        for i in range(self.n):
            for j in range(i+1, self.n):
                if self.C[i] == self.C[j]:
                    if self.P[i,j] > Wmax:
                        Wmax = self.P[i,j]
                else:
                    if self.P[i,j] < Bmin:
                        Bmin = self.P[i,j]
                a = C[i]
                b = C[j]
                if C[i] > C[j]:
                    a = C[j]
                    b = C[i]
                if self.P[i,j] > maxx[a,b]:
                    maxx[a,b] = self.P[i,j]
                if self.P[i,j] < minn[a,b]:
                    minn[a,b] = self.P[i,j]

        print "Dunn Vals", Bmin, Wmax
        print "maxvals", maxx
        print "minvals", minn
        return (Bmin/Wmax)

    def davies_bouldin(self):
        sd = np.zeros(self.k)
        for i in range(self.n):
            cl = self.C[i]
            diff = X[i,:]-self.mu[cl,:]
            sd[cl] += np.dot(diff,diff)
        for i in range(self.k):
            sd[i] /= self.nc[i]
            sd[i] = np.sqrt(sd[i])
            
        print "DBB", sd

        DB = np.zeros((self.k, self.k))
        for i in range(self.k):
            for j in range(i+1,self.k):
                diff = self.mu[i,:]-self.mu[j,:]
                DB[i,j] = (sd[i] + sd[j])/np.sqrt(np.dot(diff,diff))
                print i,j,diff, np.sqrt(np.dot(diff,diff)), sd[i]+sd[j] 
                DB[j,i] = DB[i,j]

        print DB
        DBavg = 0.0
        for i in range(self.k):
            DBavg += np.max(DB[i,:])
        DBavg /= self.k
        return (DBavg)


    def HubertInt(self):
        Y = np.zeros((self.n, self.n))
        mudist = np.zeros((self.k, self.k))
        for i in range(self.k):
            for j in range(i+1,self.k):
                diff = self.mu[i,:]-self.mu[j,:]
                dis = np.sqrt(np.dot(diff,diff))
                mudist[i,j] = dis
                mudist[j,i] = dis
        for i in range(self.n):
            for j in range(i+1, self.n):
                mud = mudist[self.C[i], self.C[j]]
                Y[i,j] = Y[j,i] = mud

        Gamma = 0.0

        N = self.n*(self.n-1)/2
        p = []
        y = []
        for i in range(self.n):
            for j in range(i+1, self.n):
                Gamma += self.P[i,j]*Y[i,j]
                p.append(self.P[i,j])
                y.append(Y[i,j])
        
        print "HubertTT", Gamma, N
        Gamma = Gamma/N
        p = np.array(p)
        y = np.array(y)
        p.shape = (N,1)
        y.shape = (N,1)
        mup = np.mean(p)
        muy = np.mean(y)
        zp = p - mup*np.ones((N,1))
        zy = y - muy*np.ones((N,1))
        
        Gamman = np.dot(zp.T,zy)/np.sqrt(np.dot(zp.T,zp)*np.dot(zy.T,zy))
        return (Gamma, Gamman[0,0])

######################
def get_dist(a,b):
    diff = a-b
    dist = np.sqrt(sum(diff*diff))
    return dist

def run_kmeans(idx):
    (kk, i) = Work[idx]
    XX = XXary[i]
    kmeans = KMeans(init='k-means++', k=kk, n_init=100).fit(XX)
    CK = kmeans.labels_
    nck = np.zeros(kk)
    for a in range(n):
        nck[CK[a]] += 1
    PK = dist(XX)
    PK.cluster_stats(CK, kk, nck,out=False)
    return PK

def get_KL_divergence(P,Q):
    l = len(P)
    kl = 0.0
    for i in range(l):
        if P[i] != 0 and Q[i] != 0:
            kl += P[i] * np.log2(P[i]/Q[i])
    return kl

#########MAIN##########

f = open(sys.argv[1]) #X.cluster.txt file
gapstatistic = False #turn on to generate gapstatistic
gengraph = False #make false after running once

k = 3 #number of true clusters
tmap = {'Iris-setosa': 0, 'Iris-versicolor':1, 'Iris-virginica':2}
#cmap = {1:0, 2:1, 3:2}
dX = []
dY = []
C = []
T = []
n = 0 #number of points

r = 0 #number of mined clusters
for l in f.readlines():
    a = l.strip().split(',')
    d = len(a)-2
    c = int(a[-1])
    t = a[-2].strip('"')
    if r < c: r = c
    C.append(c-1) #cluster id is one less than c (starts from 0)
    T.append(t)
    #if t not in tmap: 
    #    tmap[t] = k
    #    k += 1
    #if c not in cmap: 
    #    cmap[c] = r
    #    r += 1
    dY.append(tmap[t])
    dX.extend([float(x) for x in a[:-2]])
    n+= 1

#print n,d,len(dX)
X = np.array(dX)
X.shape = (n,d)
Y = np.array(dY)


if gengraph:
    #print empty nodes with names
    for i in range(n):
        print "\pnode(%0.3f, %0.3f){n%d}" % (X[i,0], X[i,1], i)

    #print internal edges
    print "\psset{linewidth=0.01pt,linecolor=lightgray}"
    for i in range(n):
        for j in range(i+1, n):
            if C[i] == C[j]:
                print "\\ncline[]{n%d}{n%d}" % (i, j)

    #print dots for nodes
    print "\psset{fillcolor=gray, linecolor=black}"
    ds = ["Bsquare", "Bo", "Btriangle"]
    for i in range(n):
        dstyle = ds[C[i]]
        print "\dotnode[dotstyle=%s](%0.3f, %0.3f){n%d}" % (dstyle, X[i,0], X[i,1], i)
    exit()



N = n*(n-1)/2 #number of point pairs

#compute CT via C and T arrays
CT = np.zeros((r,k))
for i in range(n):
    CT[C[i], tmap[T[i]]] += 1
print CT

#for key in sorted(CT2):
#    print key, CT2[key]

#compute TP, FN, FP, TN by looking at all point pairs

TP = 0
FN = 0
FP = 0
TN = 0
for i in range(n-1):
    for j in range(i+1, n):
        if T[i] == T[j] and C[i] == C[j]:
            TP += 1
        if T[i] == T[j] and C[i] != C[j]:
            FN += 1
        if T[i] != T[j] and C[i] == C[j]:
            FP += 1
        if T[i] != T[j] and C[i] != C[j]:
            TN += 1

print TP, FN, FP, TN, TP+FN+FP+TN, N

#compute using CT
TP = 0
FN = 0
FP = 0
TN = 0
nt = CT.sum(axis=0)
nc = CT.sum(axis=1)
print nc, nt
for i in range(r):
    FP += nc[i]*(nc[i]-1)/2
    for j in range(k):
        TP += CT[i][j]*(CT[i][j]-1)/2
FP = FP - TP
for j in range(k):
    FN += nt[j]*(nt[j]-1)/2
FN = FN - TP

TN = N-(TP+FN+FP)
print "TP, FN, FP, TN", TP, FN, FP, TN

#CT based measures
purity = 0.0
Fary = []
Centropy = 0.0
CEary = []
Hc = 0.0
MI = 0.0
for i in range(r):
    tji = np.argmax(CT[i,:])
    nji = np.max(CT[i,:])
    purity += nji
    Fi = (2.0*nji)/(np.sum(CT[i,:])+np.sum(CT[:,tji]))
    Fary.append(Fi)
    pi = nc[i]*1.0/n
    print "piii", pi * np.log2(pi)
    Hc += pi * np.log2(pi)
    cei = 0
    for j in range(k):
        pj = nt[j]*1.0/n
        nij = CT[i][j]
        pij = nij*1.0/n
        if pij > 0:
            Centropy += pij * np.log2(nij/nc[i])
            cei += nij/nc[i] * np.log2(nij/nc[i])
            MI += pij * np.log2((n*nij*1.0)/(nc[i]*nt[j]))
            print "i,j", i, j, pij * np.log2((n*nij*1.0)/(nc[i]*nt[j]))

    cei = -1.0*cei
    CEary.append(cei)        

Ht = 0.0
for j in range(k):
    pj = nt[j]*1.0/n
    Ht += pj * np.log2(pj)

Hc = -1.0*Hc
Ht = -1.0*Ht
NMI = MI/np.sqrt(Hc*Ht)
VI = Hc + Ht - 2*MI
print "MI, Hc, Ht, VI", MI, Hc, Ht, VI
Centropy = -1.0*Centropy
purity = purity/n
print "Fary", Fary
print "CEary", CEary

F = np.mean(Fary)

G = nx.Graph()
for i in range(r):
    for j in range(k):
        G.add_edge("c"+str(i), "t"+str(j), weight=CT[i][j])

#print G.nodes()
#print G.edges()
M = nx.max_weight_matching(G)

wM = 0
for i in range(r):
    ci = "c"+str(i)
    if ci in M:
        tj = M[ci]
        wM += G[ci][tj]['weight']
print "matching", wM, M
wM = wM*1.0/n
print "CT: Pur, F, Ent, wM, NMI:", purity, F, Centropy, wM, NMI

#Jaccard, Rand, FM
Jaccard = TP*1.0/(TP+FN+FP)
Rand = (TP+TN)/(N*1.0)
FM = np.sqrt(TP**2/((TP+FN)*(TP+FP)*1.0))
print "J, R, FM", Jaccard, Rand, FM

#verify that FM is also a dot product
bC = np.zeros((n,n))
bT = np.zeros((n,n))
for i in range(n):
    for j in range(n):
        if i == j: continue
        if C[i] == C[j]: bC[i,j]=1
        if T[i] == T[j]: bT[i,j]=1
FM2 = np.sum(bC*bT)/np.sqrt(np.sum(bC*bC) * np.sum(bT*bT))
print "FM dot", FM2

#Hubert's Gamma
bC = np.zeros((n,n))
bT = np.zeros((n,n))
for i in range(n-1):
    for j in range(i+1,n):
        if T[i] == T[j]: bT[i][j] = 1
        if C[i] == C[j]: bC[i][j] = 1

Gamma = np.sum(bT*bC)/N
print "Gamma", Gamma, TP*1.0/N

#normalized Gamma
muT = (TP + FN)/(N*1.0)
muC = (TP + FP)/(N*1.0)

Gamman = (TP/(N*1.0) - muT*muC)/np.sqrt(muT*muC*(1-muT)*(1-muC))
print "Gamman", Gamman



##############INTERNAL MEASURES
k = r #make sure that we now use the number of mined clusters as k
P = dist(X)
P.cluster_stats(C,k,nc)
print "cluster stats"
print "W,B", P.W, P.B
print "NW, NB", P.NW, P.NB
print "Wi", P.Wi
print "Bi", P.Bi
print "voli", P.voli
print "SW", P.SW
print "SB", P.SB


#BetaCV, NC, Q, Cindex
BetaCV = (P.NB * P.W)/(P.NW * P.B)
NC = 0.0
Q = 0.0
Wvv = 0.0
for i in range(k):
    Wvv += P.voli[i] #W(V,V)

for i in range(k):
    Q += (P.Wi[i]/Wvv - (P.voli[i]/Wvv)**2)
    print "QQ", P.Wi[i]/Wvv - (P.voli[i]/Wvv)**2
    NC += P.Bi[i]/P.voli[i]

(Wmin, Wmax) = P.sum_min_max(P.NW)
print "Wmin, Wmax", Wmin, Wmax
Cindex = (P.W/2.0 - Wmin)/(Wmax - Wmin)

print "BetaCV, NC, Q, Cindex", BetaCV, NC, Q, Cindex

#Silhouette Coefficient
SC = sk.metrics.silhouette_score(X, np.array(C), metric='euclidean')
print "SC", SC

SCdict = {}
SC = 0.0
for i in range(n):
    sil = P.silhouette(i,out=True)
    if C[i] not in SCdict: SCdict[C[i]] = []
    SCdict[C[i]].append(sil)
    SC += sil

SC /= n

fsc = open("/tmp/tmp_sil.txt", "w")
i = 0
for ci in sorted(SCdict):
    print "ClusterSC", ci, np.mean(SCdict[ci]), len(SCdict[ci])
    for val in sorted(SCdict[ci], reverse=True):
        print >>fsc, i, val, ci
        i += 1


Dunn = P.dunn()
DB = P.davies_bouldin()
print "SC, Dunn, DB", SC, Dunn, DB

#Huberts
Gamma, Gamman = P.HubertInt()
print "Gamma, Gamman", Gamma, Gamman

maxk = 10
numprocs = 12
repeats = 500

#2d histogram for tendency & 1d histo of the edge lenghts
bins = 5 #for 2d
bins1d = 25 #for 1d
H, xedges, yedges = np.histogram2d(X[:,0], X[:,1], bins=bins)
P = H/np.sum(H)
print "histogram X"
print P, H
DD = spatial.distance.pdist(X)
H1, xedges1 = np.histogram(DD.flatten(), bins1d)
P1 = (1.0*H1)/np.sum(H1)
print "histogram1D", P1, H1

(n,d) = X.shape
mins = np.min(X, axis=0)
maxs = np.max(X, axis=0)

KLvalsU = np.zeros(repeats)
KLvalsU1 = np.zeros(repeats)
#Puexample = np.zeros(bins)
#Pu1example = np.zeros(bins1d)
for i in range(repeats):
    #sample uniformly in data space
    XX = np.zeros((n,d))
    for a in range(d):
        XX[:,a] = np.random.uniform(mins[a], maxs[a], n)
    Hu, xe, ye = np.histogram2d(XX[:,0], XX[:,1], bins=bins)
    Pu = Hu/np.sum(Hu)
    KLvalsU[i] = get_KL_divergence(P.flatten(), Pu.flatten())
    DD = spatial.distance.pdist(XX)
    H1u, xe = np.histogram(DD.flatten(), bins1d)
    P1u = (1.0*H1u)/np.sum(H1u)
    KLvalsU1[i] = get_KL_divergence(P1, P1u)
    if i == 0:
        Puexample = Pu
        Pu1example = P1u

#print "KLvalsU", KLvalsU
#print "KLvalsU1", KLvalsU1
klbins = 10
print "Uniform", np.mean(KLvalsU), np.std(KLvalsU)
Hu, xe = np.histogram(KLvalsU, klbins)
Hu = (1.0*Hu)/np.sum(Hu)
midpts = np.diff(xe)/2.0+xe[:-1]
for i in range(klbins):
    print i, midpts[i], Hu[i]

print "Uniform1D", np.mean(KLvalsU1), np.std(KLvalsU1)
Hu, xe = np.histogram(KLvalsU1, klbins)
Hu = (1.0*Hu)/np.sum(Hu)
midpts = np.diff(xe)/2.0+xe[:-1]
for i in range(klbins):
    print i, midpts[i], Hu[i]


#2D spatial histograms
print "2D histogram"
print "edgesXY", xedges, yedges
midptsX = np.diff(xedges)/2.0+xedges[:-1]
midptsY = np.diff(yedges)/2.0+yedges[:-1]
for i in range(bins):
    for j in range(bins):
        print i*bins+j, P[i,j], Puexample[i,j]

#1D edge length histogram
print "1D histogram"
midpts = np.diff(xedges1)/2.0+xedges1[:-1]
for i in range(bins1d):
    print i, midpts[i], P1[i], Pu1example[i]


#Hopkins Statistic
sratio = 0.2
m = int(np.ceil(sratio*n))
KDT = spatial.cKDTree(X)
Hcnt = {}
maxcnt = 0
for i in range(n):
    key = str(X[i,:])
    if key not in Hcnt: Hcnt[key] = 0
    Hcnt[key] += 1
    if Hcnt[key] > maxcnt:
        maxcnt = Hcnt[key]
print "maxcnt", maxcnt

HSvals = np.zeros(repeats)
for i in range(repeats):
    #sample uniformly in data space
    RR = np.zeros((m,d))
    for a in range(d):
        RR[:,a] = np.random.uniform(mins[a], maxs[a], m)
    #select m points from X 
    SS = np.zeros((m,d))
    idx = np.random.randint(0,n,m)
    for b,j in enumerate(idx):
        SS[b,:] = X[j,:] 
    Rdist = 0.0
    Sdist = 0.0
    for j in range(m):
        rd, ri = KDT.query(RR[j,:], k=maxcnt+1)
        sd, si = KDT.query(SS[j,:], k=maxcnt+1)
        Rd = 0
        Sd = 0
        for a in range(maxcnt+1):
            if rd[a] > 0: 
                Rd = rd[a]
                break
        for a in range(maxcnt+1):
            if sd[a] > 0: 
                Sd = sd[a]
                break
        Rdist += Rd**d
        Sdist += Sd**d

    hs = Rdist/(Rdist + Sdist)
    HSvals[i] = hs

print "Hopkins", np.mean(HSvals), np.std(HSvals)
hsbins=25
HS, xe = np.histogram(HSvals, bins=hsbins)
HS = (1.0*HS)/np.sum(HS)
midpts = np.diff(xe)/2.0+xe[:-1]
for i in range(hsbins):
    print i, midpts[i], HS[i]




gapstatistic = False #turn on when needed
CHindex = False
if gapstatistic or CHindex:
    # good kmeans clustering for diff values of k on X
    CHvals = np.zeros(maxk)
    SCvals = np.zeros(maxk)
    Wk = np.zeros(maxk)
    for kk in range(1,maxk):
            kmeans = KMeans(init='k-means++', k=kk, n_init=100).fit(X)
            CK = kmeans.labels_
            nck = np.zeros(kk)
            for a in range(n):
                nck[CK[a]] += 1
            PK = dist(X)
            PK.cluster_stats(CK, kk, nck, out=False)
            Wk[kk] = PK.W
            if kk > 1:
                silvals = map(PK.silhouette, range(n))
                SCval = np.mean(silvals)
                CH = (np.trace(PK.SB)/(kk-1))/(np.trace(PK.SW)/(n-kk))
                SCvals[kk] = SCval
                CHvals[kk] = CH
                #print np.trace(PK.SB), PK.SB
                #print np.trace(PK.SW), PK.SW

    for kk in range(1,maxk):
        print "SC, CH", kk, SCvals[kk], CHvals[kk]

    for kk in range(2,maxk-2):
        print kk+1, (CHvals[kk+2]-CHvals[kk+1])-(CHvals[kk+1]-CHvals[kk]),\
        (Wk[kk]-Wk[kk+1]) - (Wk[kk+1]-Wk[kk+2])

    print Wk

#gap statistic distribution
if gapstatistic:
    (n,d) = X.shape
    mins = np.min(X, axis=0)
    maxs = np.max(X, axis=0)

    XXary = []
    Work = []
    for i in range(repeats):
        #sample uniformly in data space
        XX = np.zeros((n,d))
        for a in range(d):
            XX[:,a] = np.random.uniform(mins[a], maxs[a], n)
        if i == 0: #save uniform data once for plotting
            np.savetxt("/tmp/gapdata.txt", XX)
        XXary.append(XX)
        for kk in range(1,maxk):
            Work.append((kk,i))


    par = mp.Pool(numprocs)
    runs = par.map(run_kmeans, range(len(Work)))

    Wkb = {}
    for idx in range(len(Work)):
        kk, i = Work[idx]
        if kk not in Wkb: Wkb[kk] = []
        Wval = runs[idx].W
        Wkb[kk].append(Wval)


    for kk in Wkb:
        Wkbm = np.mean(np.log2(Wkb[kk]))
        sdk = 0.0
        gapk = (Wkbm - np.log2(Wk[kk]))
        for i in range(repeats):
            sdk += (np.log2(Wkb[kk][i]) - Wkbm)**2
        sdk /= repeats
        sdk = np.sqrt(sdk)
        print "GAP stats", kk, gapk, sdk, gapk-sdk, np.log2(Wk[kk]), Wkbm, Wk[kk]



