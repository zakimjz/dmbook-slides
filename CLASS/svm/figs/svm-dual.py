#!/usr/bin/env python

import sys
import random
import numpy as np

def phi2h(x): #homogeneous quadratic mapping
    r = np.zeros(d*(d+1)/2)
    k = 0
    for i in range(d):
        for j in range(i,d):
            if i==j: r[k] = x[i]*x[j]
            else: r[k] = np.sqrt(2)*x[i]*x[j]
            k += 1
    return r

kernel = "linear" #default linear kernel
loss = "quadratic" #default, choose quadratic or hinge loss

fname = sys.argv[1]
C = np.float(sys.argv[2])
if len(sys.argv) > 3:
    kernel = sys.argv[3]
if len(sys.argv) > 4:
    loss = sys.argv[4]

print "params", fname, C, kernel, loss

D = np.loadtxt(fname)
(n, d) = np.shape(D) #get input dimensions

X = D[:,0:2]
X = np.append(X, np.ones((n,1)), 1) #add extra column of ones
Y = D[:,2] #classes

if kernel == "linear":
    K = np.dot(X, X.T)
    if loss == "quadratic":
        delta = np.diag(np.tile(1.0/(2*C), n))
        K = K+delta
elif kernel == "quadratic":
    K = np.dot(X, X.T)
    K = K*K #square each element
    if loss == "quadratic":
        delta = np.diag(np.tile(1.0/(2*C), n))
        K = K+delta

eta = 1/np.diag(K)
alpha = np.zeros(n)

#start iterative stochastic gradient ascent
eps = 0.0001
t = 0
err = 1
while err > eps:
    alpha_prev = alpha.copy()
    idx = range(n)
    #random.shuffle(idx)
    for k in idx:
        ay = alpha*Y
        ayK = np.dot(ay, K[:,k])
        alpha[k] += eta[k] * (1 - Y[k]*ayK)
        if alpha[k] < 0: alpha[k] = 0
        elif alpha[k] > C: alpha[k] = C
        
    t += 1
    err = np.linalg.norm(alpha-alpha_prev)
    print t, err

print "The support vectors are:"
for (i, ai) in enumerate(alpha):
    if ai > 0: print i, ai
print "number of support vectors:", len(alpha[alpha > 0])

if kernel == "linear":
    ay = alpha*Y
    w = np.sum(X.T*ay, axis=1)
    print "hyperplane:", w
else: #quadratic homogeneous
    Z = np.array([phi2h(x) for x in X]) #row-wise mapping to phi
    ay = alpha*Y
    w = np.sum(Z.T*ay, axis=1)
    print "w:", w
    z = phi2h(np.ones(d))
    print "h:", w*z


