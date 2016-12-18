#!/usr/bin/env python

import networkx as nx

f = open("iris-layout.txt", "rU")

G = nx.Graph()
for l in f.readlines():
    a = l.strip().split()
    if a[0] == "e":
        G.add_edge(a[1], a[2])

print nx.average_clustering(G)
