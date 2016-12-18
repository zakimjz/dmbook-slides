#!/usr/bin/env python

import sys

dfile = sys.argv[1]
cols = [int(x) for x in sys.argv[2].strip().split()]

f = open(dfile, "rU")

for l in f.readlines():
    a = l.strip().split()
    for c in cols:
        print a[c],
    print
