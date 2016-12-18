#!/usr/bin/env python

from numpy import *
set_printoptions(suppress=True)

f = open("iris-slw.dat", "rU")
of = open("iris-slw-c.dat", "w")

D = array([map(float, l.split()) for l in f])
d= 2
n= 150
D.shape = n,d
m = mean(D,axis=0)
print m

#center the data
Z = D-m
savetxt(of,Z,fmt="%0.4f")

x = argmax([vdot(Z[i,:], Z[i,:]) for i in range(n)])
r = vdot(Z[x,:], Z[x,:])
print sqrt(r), Z[x,:]

l = Z.min()
u = Z.max()
print "range", l, u
