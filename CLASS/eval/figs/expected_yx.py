#!/usr/bin/env python

f = open('iris-PC.txt')

H = {}
Y = {}
for l in f.readlines():
    a = l.strip().split()
    key = "%0.2f %0.2f" % (float(a[0]), float(a[1]))
    if key in H:
        H[key] += 1
        Y[key].append(a[2])
    else:
        H[key] = 1
        Y[key] = [a[2]]

for k in H:
    if H[k] > 1:
        print k, H[k], Y[k]
