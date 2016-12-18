#!/usr/bin/env python

import numpy as np

u1 = np.matrix([-0.390, 0.089, -0.916])
x = np.matrix([-0.343, -0.754, 0.241])
a = x*u1.T
x2 = a*u1
print x2
e = x-x2
print e
u2 = np.matrix([-0.639, -0.742, 0.2])
u3 = np.matrix([-0.663, 0.664, 0.346])
e2 = 0.828*u2 - 0.19*u3
print e2

#projection matrix
P1  =u1.T*u1
print P1
xp = P1*x.T
print xp.T*e.T

print u2.T*u2
P2 = u1.T*u1 + u2.T*u2
print P2


U2 = np.bmat([u1.T,u2.T])
print "U2", U2
print U2*U2.T

