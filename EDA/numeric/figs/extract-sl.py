#!/usr/bin/env python

f = open('iris.txt', 'rU')

cnt = {}
for l in f.readlines():
    a = l.strip().split(',')
    if a[0] not in cnt: cnt[a[0]] = 0
    cnt[a[0]] += 0.05
    print a[0], cnt[a[0]]
