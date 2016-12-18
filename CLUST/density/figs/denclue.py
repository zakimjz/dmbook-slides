#!/usr/bin/env python

import sys, os
import random
import numpy as np
import networkx as nx
import multiprocessing as mp
from scipy import spatial
from time import time

#helper routines
clmap = {"Iris-versicolor": 0, "Iris-setosa": 1, 
         "Iris-virginica": 2}
def class2int(c):
    if fname == "iris.txt":
        c = c.strip('"')
        return clmap[c]
    elif fname == "iris-PC.txt":
        return int(c)

#std mv norm density
def dmvnorm_std(x):
    part1 = np.power(2*np.pi, d/2)
    prod = np.dot(x, x)
    part2 = np.exp(-0.5*prod)
    res =  part2/part1
    return res



#update denclue2.0 rule
def update_x(x):
	if use_nn == 1:
		radius = 5*h
		Neps = kdtree.query_ball_point(x, radius)
		Z = (x-X[Neps,])/h
		Kvals = map(dmvnorm_std, Z)
		xnew = np.sum(Kvals*X[Neps,].T, axis=1)/np.sum(Kvals)
	else:
		Z = (x-X)/h
		Kvals = map(dmvnorm_std, Z)
		xnew = np.sum(Kvals*X.T, axis=1)/np.sum(Kvals)

	return xnew


#update denclue1.0 rule (gradient based)
def gradient_update_x(x):
    step = 0.05
    Z = (x-X)/h
    Kvals = map(dmvnorm_std, Z)
    Kvals = [val*h*-1.0 for val in Kvals]
    grad = np.sum(Kvals*Z.T, axis=1)
    grad = grad/(n*h**(d+2))
    xnew = x + step*grad
    return xnew

def hessian_update_x(x):
    Z = (x-X)/h
    Kvals = map(dmvnorm_std, Z)
    Hi = np.eye(d,d)*np.sum(Kvals)/(n*h**(d+2))
    H = np.dot(Kvals*Z.T, Z)/(n*h**(d+4))
    H = H+Hi
    
    Kvals = [val*h*-1.0 for val in Kvals]
    grad = np.sum(Kvals*Z.T, axis=1)
    grad = grad/(n*h**(d+2))

    xnew = x + np.dot(np.linalg.inv(H), grad)
    return xnew



def find_attractor(idx):
        notconverged = True
        t = 0
        x = X[idx,:]
        while(notconverged):
            if denclue == "2.0": 
                xnew = update_x(x)
            else: 
                xnew = gradient_update_x(x)
                #xnew = hessian_update_x(x)

            diff = np.linalg.norm(xnew-x)
            notconverged = diff > eps
            t = t+1
            x = xnew

        print idx, t, x, os.getpid()
        return x

############main#############

h = 0.5 #spread
numprocs = 8 #how many processes to use
denclue = "2.0" #use denclue2.0 update rule
use_nn = 1
#run as: denclue fname classattr-pos h p version(2.0/1.0)
fname = sys.argv[1]
classattr = int(sys.argv[2])
if len(sys.argv) > 3:
    h = float(sys.argv[3])
if len(sys.argv) > 4:
    numprocs = int(sys.argv[4])
if len(sys.argv) > 5:
    denclue = sys.argv[5]
if len(sys.argv) > 6:
    use_nn = int(sys.argv[6]) # 0 or 1


D = np.loadtxt(fname,delimiter=",", 
               usecols = (0,1,classattr),
               converters={classattr:class2int})
(n, d) = np.shape(D) #get input dimensions
#d = d-1
d=2
X = D[:,0:d]
Y = D[:,d]

#compute only over nearest neighbors
kdtree = spatial.KDTree(X)

start = time()
eps = 0.001 #convergence threshold
#run density finding in parallel
par = mp.Pool(numprocs)
attractors = par.map(find_attractor, range(n))
elapsed = time()-start
print "elapsed time", elapsed

##sequential
#attractors = []
#for i in range(n):
    #at = find_attractor(i)
    #attractors.append(at)

#find unique attractors, assign points to clusters
#use two digits after decimal for proximity
clusters = {}
for i in range(n):
    hashL = ["%0.2f"%x for x in attractors[i]]
    hashval = " ".join(hashL)
    if hashval not in clusters: clusters[hashval] = []
    clusters[hashval].append(i)

#merge attractors/clusters that are close
#use connected components method
ceps = 0.1 #merge two clusters if attractors are closer than ceps
Uattrs = [hval for hval in clusters]
G = nx.Graph()
for i in range(len(Uattrs)):
    G.add_node(i)
    for j in range(i+1,len(Uattrs)):
        ai = np.array([float(x) for x in Uattrs[i].split()])
        aj = np.array([float(x) for x in Uattrs[j].split()])
        if np.linalg.norm(ai-aj) <= ceps:
            G.add_edge(i,j)

print G.size(), G.order()

CC = nx.connected_components(G)
cid = 0 #cluster id
for comp in CC:
    print "cluster", cid, comp
    size = 0
    pts = []
    mean = np.zeros(d)
    for v in comp:
        hval = Uattrs[v]
        print "\t", hval, len(clusters[hval])
        size += len(clusters[hval])
        pts.extend(clusters[hval])
        mui = np.array([float(x) for x in hval.split()])
        mean += mui
    mean = mean/len(comp)
    print "\tmean", mean
    print "\tpts", sorted(pts)
    print "\tsize", size
    cid += 1
