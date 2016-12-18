#!/usr/bin/env python
import numpy as np

#helper routines
clmap = {"C0": 0, "C1": 1, "C2": 2, "C3":3}
def class2int(c):
    c = c.strip('"')
    return clmap[c]

######MAIN#######
o = open("test1.dat", "w")
oa = open("test1a.dat", "w")
ob = open("test1b.dat", "w")
oc = open("test1c.dat", "w")

X = np.loadtxt("test1.elki.txt", converters={3:class2int})
(n,d) = np.shape(X)

mu = np.mean(X, axis=0)
print "mean", mu
mins = map(np.min, [X[:,0:d-1]])
maxs = map(np.max, [X[:,0:d-1]])
print mins
print maxs

for c in range(3):
    Z = X[X[:,-1] == c, :]
    mu = np.mean(Z, axis=0)
    print "mean before", c, mu
    minZ = map(np.min, [Z[:,0:d-1]])
    maxZ = map(np.max, [Z[:,0:d-1]])
    print minZ
    print maxZ


scale = 1
for i in range(n):
    if i < 100: oo = oa
    elif i < 200: oo = ob
    else: oo = oc
    for j in range(d-1):
        #X[i,j] = (X[i,j] - mins[j])*scale/(maxs[j]-mins[j])
        X[i,j] = X[i,j]*scale
    
    print >>o, " ".join(["%0.3f"%x for x in X[i,0:d-1]])
    print >>oo, " ".join(["%0.3f"%x for x in X[i,0:d-1]])

mins = map(np.min, [X[:,0:d-1]])
maxs = map(np.max, [X[:,0:d-1]])
print mins
print maxs

for c in range(3):
    Z = X[X[:,-1] == c, :]
    mu = np.mean(Z, axis=0)
    print "mean", c, mu
    minZ = map(np.min, [Z[:,0:d-1]])
    maxZ = map(np.max, [Z[:,0:d-1]])
    print minZ
    print maxZ



