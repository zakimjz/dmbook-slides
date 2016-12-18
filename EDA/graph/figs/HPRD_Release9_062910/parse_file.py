#!/usr/bin/env python

f = open("BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt", "rU")

edges = {}
for l in f.readlines():
    a = l.strip().split()
    p1 = int(a[1])
    p2 = int(a[4])
    if p2 < p1: (p1,p2) = (p2,p1)
    if (p1 != p2): 
        if (p1,p2) not in edges:
            edges[(p1,p2)] = 1

for (p1,p2) in sorted(edges):
    print p1, p2

fo = open("hprd.dot", "w")
print >>fo, "graph G {"
for (p1,p2) in sorted(edges):
    print >>fo, p1, "--", p2
print >>fo, "}"
