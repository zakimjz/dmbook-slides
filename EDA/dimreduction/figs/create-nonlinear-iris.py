#!/usr/bin/env python

import numpy as np
import numpy.linalg as la

X = np.loadtxt('iris-slwpl-c.dat')

V1 = 0.2*X[:,0]**2 + X[:,1]**2 + 0.1*X[:,0]*X[:,1]
V2 = X[:,1]
Y = np.transpose(np.vstack((V1,V2)))
m = np.mean(Y, axis=0)
Z = Y-m

print np.min(Z[:,0]), np.max(Z[:,0])
print np.min(Z[:,1]), np.max(Z[:,1])
np.savetxt('iris-2d-nonlinear.dat', Z, fmt="%0.4f")

nr, nc = np.shape(Z)
S = np.cov(np.transpose(Z))*(nr-1)/nr
print S

w, v = la.eigh(S)
print "eigvals, eigvecs"
print w
print v

U = np.asmatrix([v[:,1],v[:,0]]).T
print U
A = np.dot(U.T,Z.T).T
#print A
np.savetxt('iris-nonlinear-PCA.dat', A, fmt="%0.4f")
print np.min(A[:,0]), np.max(A[:,0])
print np.min(A[:,1]), np.max(A[:,1])

