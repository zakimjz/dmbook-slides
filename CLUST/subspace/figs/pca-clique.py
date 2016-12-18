#!/usr/bin/env python
import numpy as np

def pca(X):
    (n,d) = np.shape(X)
    mu = np.sum(X,axis=0)
    mu = mu/n
    Z = X - mu
    S = np.dot(Z.T,Z)/n 
    L,U = np.linalg.eig(S)
    #sort according to eigvals, since eig() does not do that
    iv = enumerate(L)
    siv = sorted(iv, key=lambda x: x[1], reverse=True)
    idx  = [i for (i,v) in siv]
    L= L[idx]
    U= U[idx,:]
    return (L,U)

def project(X, U):
    (n,d) = np.shape(X)
    mu = np.sum(X,axis=0)
    mu = mu/n
    Z = X - mu
    PZ = np.dot(U.T, Z.T)
    PZ = PZ.T
    mins = map(np.min, [PZ[:,0], PZ[:,1], PZ[:,2]])
    maxs = map(np.max, [PZ[:,0], PZ[:,1], PZ[:,2]])
    abss = map(np.max, np.abs([PZ[:,0], PZ[:,1], PZ[:,2]]))
    print "mins", mins
    print "maxs", maxs
    Ua = U*abss
    A = np.diagflat(abss)
    print np.dot(U,A)
    print "uam", np.dot(U,A).T+mu
    print "uam2", -1*np.dot(U,A).T+mu
    dd = np.sum(Ua, axis=1)
    print dd
    print "mean", mu
    print "U", U
    print "cov", np.cov(X.T, bias=1)
    print "mu+dd", dd+mu
    print "vecs", (-2*Ua).T

    return PZ


####MAIN#######
X = np.loadtxt("test1.dat")
(n,d) = np.shape(X)

###PCA
(L, U) = pca(X) #Z is the centered data
print "X", L
print U
PZ = project(X, U)
oo = open("test1.pca.dat", "w")
#print >>oo, "%% %0.3f %0.3f %0.3f"% (L[0], L[1], L[2])
for i in range(n):
    print >>oo, "%0.3f %0.3f %0.3f"% (PZ[i,0], PZ[i,1], PZ[i,2])

Xa = np.loadtxt("test1a.dat")
Xb = np.loadtxt("test1b.dat")
Xc = np.loadtxt("test1c.dat")

(L, U) = pca(Xa)
print "Xa", L
print np.cumsum(L[::-1])/np.sum(L)
print U
PZ = project(Xa, U)


(L, U) = pca(Xb)
print "Xb", L
print np.cumsum(L[::-1])/np.sum(L)
print U
PZ = project(Xb, U)

(L, U) = pca(Xc)
print "Xc", L
print np.cumsum(L[::-1])/np.sum(L)
print U
PZ = project(Xc, U)

###CLIQUE
w = 1
#clique(X)


