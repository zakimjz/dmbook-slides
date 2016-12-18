#!/usr/bin/env python

import networkx as nx
f = open("karate.gml")

G=nx.Graph()
for l in f.readlines():
    a = l.strip().split()
    if len(a) < 1: next
    if a[0] == "source":
        s = a[1]
    elif a[0] == "target":
        t = a[1]
        if s > t: (s,t) = (t,s)
        G.add_edge(s,t)

A = nx.adj_matrix(G)
n = G.order()
for i in range(0,n):
    for j in range(0,n):
        print A[i,j],
    print
