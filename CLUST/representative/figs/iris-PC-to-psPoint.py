#!/usr/bin/env python

import sys

printpspoints = True
f = open(sys.argv[1])
printpspoints = int(sys.argv[2])
foutname = sys.argv[3]

clusters = {}
purity = {}
pnum = 1
for l in f.readlines():
    a = l.strip().split(",")
    (x,y) = (float(a[0]), float(a[1]))
    tid = a[2] #true class id
    if len(a) > 3: cid = int(a[3]) #predicted cluster id
    else: cid = 1

    if cid not in purity: purity[cid] = {}
    if tid not in purity[cid]: purity[cid][tid] = 0
    purity[cid][tid] += 1
    if cid not in clusters: clusters[cid] = []
    clusters[cid].append((pnum, x, y, tid, cid))
    pnum += 1

#get correct max class
classid = {}
for cid in purity:
    for tid in sorted(purity[cid], key=purity[cid].get, reverse=True):
        classid[cid] = tid
        break

print classid

fo = {}
fw = {}
for k in clusters:
    fo[k] = open(foutname+"-C"+str(k)+".tex", "w")
    fw[k] = open(foutname+"-W"+str(k)+".tex", "w")

for k in clusters:
    for (p,x,y, tid, cid) in clusters[k]:
        fname = fo[k]
        if tid != classid[cid]:
            fname = fw[k]
        if printpspoints:
            print >>fname, "\psPoint(%3.2f,%3.2f,0){p%d}"%(x,y,p)
            print >>fname, "\psdots[](p%d)" % (p)
        else:
            print >>fname, "\psdot[](%3.2f,%3.2f)" % (x,y)

#print(purity)
