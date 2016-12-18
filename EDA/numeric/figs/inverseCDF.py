#!/usr/bin/env python

f = open('iris.txt', 'rU')

sl = {}
n=0.0
for l in f.readlines():
    a = l.strip().split(',')
    if a[0] not in sl: sl[a[0]]=0
    sl[a[0]] += 1
    n+=1

cumsum = 0
for k in sorted(sl):
    cumsum += sl[k]
    sl[k] = cumsum

cdf = open('sl-CDF.dat', 'w')
icdf = open('sl-iCDF.dat', 'w')

prevP = -1
for k in sorted(sl):
    if prevP != -1:
        print >> icdf, prevP, k
    print >>icdf, sl[k]/n, k
    if prevP != -1: 
        print >>cdf, k, prevP 
    print >>cdf, k, sl[k]/n
    prevP = sl[k]/n

