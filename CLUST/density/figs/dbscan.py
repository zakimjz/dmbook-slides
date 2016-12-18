#!/usr/bin/env python
import sys, os
import numpy as np
from scipy import spatial
import multiprocessing as mp
import matplotlib.pyplot as plt

#helper routines
clmap = {"Iris-versicolor": 0, "Iris-setosa": 1, 
         "Iris-virginica": 2}
def class2int(c):
    if fname == "iris.txt":
        c = c.strip('"')
        return clmap[c]
    elif fname == "iris-PC.txt":
        return int(c)
    else:
        return 0


def mark_reachable(pt, cluster, Marked, Cores, Borders, Neps):
    #print "Neighbors", pt, Neps[pt]
    for i in Neps[pt]:
        cluster[i] = 1  #add point
        if i in Cores: #i is core point
            if i not in Marked:
                Marked[i] = True
                mark_reachable(i, cluster, Marked,\
                               Cores, Borders, Neps)
        else:
            Borders[i] = 1


def classify_pts(X, Neps):
    #compute core points
    Cores = {}
    for i in range(n):
        if len(Neps[i]) >= minpts:
            Cores[i] = 1

    cid = 0
    Clusters = {}
    Borders = {}
    Marked = {}
    #find the clusters
    for core in Cores:
        if core not in Marked:
            Marked[core] = True
            Clusters[cid] = {} #create new cluster with core
            Clusters[cid][core] = 1
            mark_reachable(core, Clusters[cid], Marked,\
                           Cores, Borders, Neps)
            cid += 1

    return (Clusters, Cores, Borders)


def getNeps(i):
    res = kdtree.query_ball_point(X[i,:], eps)
    return res


######### main ##########
#run as: dbscan fname classattr-pos eps minpts procs

fname = sys.argv[1] #only first two attrs used
classattr = int(sys.argv[2])
eps = float(sys.argv[3])
minpts = int(sys.argv[4])

numprocs = 8 #how many processes to use
if len(sys.argv) > 5:
    numprocs = int(sys.argv[5])

D = np.loadtxt(fname,delimiter=",", 
               usecols = (0,1,classattr),
               converters={classattr:class2int})
(n, d) = np.shape(D) #get input dimensions
d=2
X = D[:,0:d]

#compute Neps and find core points
kdtree = spatial.KDTree(X)
Neps = []
par = mp.Pool(numprocs)
Neps = par.map(getNeps, range(n))
#for i in range(n):
    #res = getNeps(i)
    #Neps.append(res)

sys.setrecursionlimit(n)
(clusters, cores, borders) = classify_pts(X, Neps)

totsize = 0
noise = {x:1 for x in range(n)}
colors = 'bgrcmyk'
symbols = 'o*+xsD1234'
ofname = open(fname+".cluster.txt", "w")
print >>ofname, "%", fname, "minpts", minpts, "eps", eps
cnts = 0
Ccnts = 0
Bcnts = 0
for k in sorted(clusters):
    print "C", k, ":",
    pts = clusters[k].keys()
    print len(pts)#, pts
    Cpts = [x for x in pts if x in cores]
    cnts += len(Cpts)
    Ccnts = cnts
    print "cnts ", cnts
    Bpts = [x for x in pts if x not in cores]
    cnts += len(Bpts)
    Bcnts = cnts
    print "cnts ", cnts

    print >>ofname, "%", k, len(Cpts), len(Bpts), Ccnts, Bcnts
    for x in Cpts:
        strx = "%0.2f, %0.2f" %( X[x,0], X[x,1])
        print >>ofname, strx
    for x in Bpts:
        strx = "%0.2f, %0.2f" %( X[x,0], X[x,1])
        print >>ofname, strx
    plt.plot(X[pts,0], X[pts,1],\
             colors[k%len(colors)]+symbols[k%len(symbols)])
    for v in pts: 
        if v in noise: del(noise[v])
    totsize += len(pts)

noisep = noise.keys()
print >>ofname, "% noise", len(noisep)
for x in noisep:
    strx = "%0.2f, %0.2f" %( X[x,0], X[x,1])
    print >>ofname, strx

plt.plot(X[noisep,0], X[noisep,1], 'c.')
plt.show()
print totsize
print "cores", len(cores)#, cores.keys()
print "borders", len(borders)#, borders.keys()
print "noise", len(noise)#, noise.keys()

