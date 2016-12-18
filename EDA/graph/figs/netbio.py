#!/usr/bin/env python

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from math import log
import operator

f = open("hprd.txt", "rU")
fd = open("hprd-deg.dat", "w")
fc = open("hprd-ck.dat", "w")

#PPI graph
G = nx.Graph()
for l in f.readlines():
    a = l.strip().split()
    if a[0] != a[1]:
        G.add_edge(a[0],a[1])
print G.size(), G.order()

#nx.draw_spring(G)
#plt.savefig("G.pdf")

#diameter
#splen = nx.all_pairs_shortest_path_length(G)
#diameter = -1
#for n in G.nodes():
    #for n2 in G.nodes():
        #if n in splen and n2 in splen[n] and splen[n][n2] > diameter:
            #diameter = splen[n][n2]
#print "diameter =", diameter

#compute degree distribution
Pk = {} #deg dist
for n in G.nodes():
    dn = G.degree(n)
    if dn not in Pk: Pk[dn] = 0
    Pk[dn] += 1

D = [G.degree(n) for n in G.nodes()]
print "max deg node", max(D), D.index(max(D))

#print Pk
R = []
rsum = 1
for k in sorted(Pk):
    R.append(log(rsum,2))
    print >>fd, k, 1.0*Pk[k]/G.order(), rsum
    rsum -= 1.0*Pk[k]/G.order()


xs = 3
xe = 270 #max degree to consider for fitting
#md = max(Pk)
#print md

K = [log(k,2) for k in sorted(Pk)]
P = [log(1.0*Pk[k]/G.order(),2) for k in sorted(Pk)]
#print K
#print P
plt.plot(K,P,'ro')
plt.ylabel('log P(k): probability')
plt.xlabel('log k: degree')
m = np.polyfit(K[xs:xe], P[xs:xe], 1)
#print m
pfit = np.polyval(m,K[xs:xe])
s = "slope=%0.2f"%(m[0])
plt.text(2.5, -4, s)
plt.plot(K[xs:xe],pfit,'b-')
plt.savefig("degdist.pdf")
plt.close()

plt.plot(K,R,'ro')
plt.ylabel('log P(k): probability')
plt.xlabel('log k: degree')
m = np.polyfit(K[xs:xe], R[xs:xe], 1)
#print m
pfit = np.polyval(m,K[xs:xe])
s = "slope=%0.2f"%(m[0])
plt.text(2.5, -4, s)
plt.plot(K[xs:xe],pfit,'b-')
plt.savefig("degcdist.pdf")
plt.close()

#for i in range(xs,xe):
    #    print i, K[i], P[i]

#clustering coefficient
Ck = {}
sumCn = 0.0
for n in G.nodes():
    dn = G.degree(n)
    Nn = G.neighbors(n)
    Ne = 0.0
    for j in Nn:
        Nj = G.neighbors(j)
        Common = filter(lambda(x): x in Nn, Nj)
        Ne += len(Common)
    if dn > 1:
        if dn not in Ck: Ck[dn] = []
        Cn = Ne/(dn*(dn-1))
        Ck[dn].append(Cn)
        if Cn > 1: print "error", n, Cn
        sumCn += Cn

print "average clustering coefficient =", sumCn/G.order(), nx.average_clustering(G)


csum = 0
for k in sorted(Ck):
    nn = len(Ck[k])
    ss = sum(Ck[k])
    csum += ss*1.0/nn

R = []
rsum = csum
for k in sorted(Ck):
    R.append(log(rsum,2))
    nn = len(Ck[k])
    ss = sum(Ck[k])
    print >>fc, k, ss*1.0/nn
    rsum -= ss*1.0/nn


K = [log(k,2) for k in sorted(Ck)]
P = [log(1.0*sum(Ck[k])/len(Ck[k]),2) for k in sorted(Ck)]
#plt.clear()
xs=3
xe=len(K)
plt.plot(K,P,'ro')
plt.ylabel('log C(k)')
plt.xlabel('log k: degree')
m = np.polyfit(K[xs:xe], P[xs:xe], 1)
#print m
pfit = np.polyval(m,K[xs:xe])
s = "slope=%0.2f"%(m[0])
plt.text(2.5, -4, s)
plt.plot(K[xs:xe],pfit,'b-')
plt.savefig("Ckdist.pdf")
plt.close()
print xe

plt.plot(K,R,'ro')
plt.ylabel('log C(k)')
plt.xlabel('log k: degree')
m = np.polyfit(K[xs:xe], R[xs:xe], 1)
#print m
pfit = np.polyval(m,K[xs:xe])
s = "slope=%0.2f"%(m[0])
plt.text(2.5, -4, s)
plt.plot(K[xs:xe],pfit,'b-')
plt.savefig("Ckcdist.pdf")
plt.close()

