#!/usr/bin/env python
import sys

f = open(sys.argv[1])

n = 0
PN = []
for l in f.readlines():
    p = l.strip().split()
    pname = "p"+str(n)
    PN.append(pname)
    print "\psPoint("+p[0]+","+p[1]+",0){"+pname+"}"
    n += 1

print "\psdots[]",
for pn in PN:
    print "("+pn+")",
print
