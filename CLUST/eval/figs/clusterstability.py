#!/usr/bin/env python

import sys
import numpy as np
import random
from sklearn.cluster import KMeans
import multiprocessing as mp
from multiprocessing import Lock

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

def get_cluster_assignment(mu, Z):
    (N,d) = Z.shape
    (k,d) = mu.shape
    C = np.zeros(N)
    nc = np.zeros(k)
    for i in range(N):
        dists = [get_dist(mu[j,:], Z[i,:]) for j in range(k)]
        C[i] = np.argmin(dists)
        nc[C[i]] += 1

    return (C, nc)

def get_VI_dist(idx):
    #print "VI", idx, WorkD[idx]
    (kk,(PKa, idxa),(PKb,idxb)) = WorkD[idx]

    #get merged multiset dataset, using max counts for each point
    cnta = np.zeros(n, dtype=int)
    cntb = np.zeros(n, dtype=int)
    for i in range(n):
        cnta[idxa[i]] += 1
        cntb[idxb[i]] += 1
    Z = []
    for i in range(n):
        maxc = min(cnta[i], cntb[i])
        for j in range(maxc):
            Z.append(X[i,:])
    Z = np.array(Z)
    (N,d) = Z.shape

    #print N, Z
    Ca, nca = get_cluster_assignment(PKa.mu, Z)
    Cb, ncb = get_cluster_assignment(PKb.mu, Z)
    #print nca, Ca
    #print ncb, Cb

    #compute VI
    Ha = 0.0
    Hb = 0.0
    for j in range(kk):
        pja = nca[j]*1.0/N
        pjb = ncb[j]*1.0/N
        if pja > 0: Ha += (pja * np.log2(pja))
        if pjb > 0: Hb += (pjb * np.log2(pjb))

    nab = np.zeros((kk,kk))
    for i in range(N):
        nab[Ca[i],Cb[i]] += 1

    #print nab
    Hab = 0.0
    for i in range(kk):
        for j in range(kk):
            pij = nab[i,j]*1.0/N 
            if pij > 0: Hab += (pij * np.log2(pij))

    VI = Ha + Hb - 2*Hab

    #compute FM
    nij2 = np.sum([nab[i,j]**2 for i in range(kk) for j in range(kk)])
    ni2 = np.sum([nca[i]**2 for i in range(kk)])
    mj2 = np.sum([ncb[i]**2 for i in range(kk)])
    TP = 0.5*(nij2 - N)
    FN = 0.5*(mj2 - nij2)
    FP = 0.5*(ni2 - nij2)
    FM = TP/np.sqrt((TP+FN)*(TP+FP))

    #silA = np.mean([PKa.silhouette(i) for i in range(n)])
    #silB = np.mean([PKb.silhouette(i) for i in range(n)])
    #sil = np.abs(silA-silB)
    #lock.acquire()
    #print idx, dist, N, d, Ha, Hb, Hab
    #lock.release()
    
    return (VI, FM)


#########MAIN##########

f = open(sys.argv[1]) #X.cluster.txt file

X = []
n = 0 #number of points

for l in f.readlines():
    a = l.strip().split(',')
    d = len(a)-2
    X.extend([float(x) for x in a[:-2]])
    n+= 1

#print n,d,len(dX)
X = np.array(X)
X.shape = (n,d)
print n,d

maxk = 11
numprocs = 12
repeats = 100

#cluster stability based on bootstrapping
XXary = []
XXidx = []
Work = []
lock = Lock()
for i in range(repeats):
    #generate a sample of size n with replacement
    ridx = [random.choice(range(n)) for j in range(n)]
    XX = np.zeros((n,d))
    for j in range(n):
        XX[j,:] = X[ridx[j]]
    XXary.append(XX)
    XXidx.append(ridx)
    for kk in range(2,maxk):
        Work.append((kk,i))

par = mp.Pool(numprocs)
print "run kmeans"
runs = par.map(run_kmeans, range(len(Work)))

Cvals = {}
for idx in range(len(Work)):
    kk, i = Work[idx]
    #print "CVAL", kk, i
    if kk not in Cvals: Cvals[kk] = []
    Cvals[kk].append((runs[idx], XXidx[i]))

print "get mean stability"
WorkD = []
for kk in Cvals:
    WorkD.extend([(kk,a,b) for (i,a) in enumerate(Cvals[kk])\
            for (j,b) in enumerate(Cvals[kk])\
            if j>i])

print "WorkD", len(WorkD)

parD = mp.Pool(numprocs)
dists = parD.map(get_VI_dist, range(len(WorkD)))
vistab = {}
fmstab = {}
for i in range(len(WorkD)):
    (kk,a,b) = WorkD[i]
    if kk not in vistab: vistab[kk] = []
    if kk not in fmstab: fmstab[kk] = []
    (vi, fm) = dists[i]
    vistab[kk].append(vi)
    fmstab[kk].append(fm)

for kk in vistab:
    print "stability: %d %0.3f %0.3f %0.3f %0.3f" %\
            (kk,\
            np.mean(vistab[kk]), np.std(vistab[kk]),\
            np.mean(fmstab[kk]), np.std(fmstab[kk]))
