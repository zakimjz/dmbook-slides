#!/usr/bin/env python

import sys
import networkx as nx
import time

def read_graphs(fname):
    #read graphs
    f = open(fname, "rU")
    eL = {}
    vL = {}
    V = {}
    N = {}
    D = {}
    for l in f.readlines():
        a = l.strip().split()
        if a[0] == 't': #graph id 
            gid = int(a[2])
            G = nx.Graph(gid=gid)
            D[gid] = G
        elif a[0] == 'v': #vertices
            vid = int(a[1])
            vlbl = a[2]
            G.add_node(vid, L=vlbl)
            if gid not in V: V[gid] = []
            V[gid].append(vid)
            vL[(gid,vid)] = vlbl
        elif a[0] == 'e': #edges
            vi = int(a[1])
            vj = int(a[2])
            elbl = 'e'
            if len(a) > 3: 
                elbl = a[3]
            G.add_edge(vi, vj, L=elbl)
            eL[(gid,vi,vj)] = elbl
            eL[(gid,vj,vi)] = elbl
            if (gid,vi) not in N: N[(gid,vi)] = []
            if (gid,vj) not in N: N[(gid,vj)] = []
            N[(gid,vi)].append(vj)
            N[(gid,vj)].append(vi)
    
    return D, vL, eL, V, N

######## main #########
fname = sys.argv[1] #graph filename

DB, VL, EL, V, N = read_graphs(fname)

rep = 1000

st = time.time()
for i in range(rep):
    for gid in DB:
        G = DB[gid]
        for u in G.nodes():
            lbl = G.node[u]['L']
            for v in G.neighbors(u):
                lbl = G.edge[u][v]['L']
et = time.time()
print "elapsed NX", et-st

st = time.time()
for i in range(rep):
    for gid in V:
        for u in V[gid]:
            lbl = VL[(gid,u)]
            #print u, gid, N[(gid,u)]
            if (gid,u) in N:
                for v in N[(gid,u)]:
                    lbl = EL[(gid,u,v)]
et = time.time()
print "elapsed", et-st
