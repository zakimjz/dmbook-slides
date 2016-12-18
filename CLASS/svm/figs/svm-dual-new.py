#!/usr/bin/env python

import sys
import random
import numpy as np

def phi2h(x): #homogeneous quadratic mapping
    d = len(x)
    r = np.zeros(d*(d-1)/2 + 1)
    k = 0
    for i in range(d-1):
        for j in range(i,d-1):
            if i==j: r[k] = x[i]*x[j]
            else: r[k] = np.sqrt(2)*x[i]*x[j]
            k += 1
    r[k] = x[d-1]*x[d-1] #last column is constant
    return r

def phi2ih(x): #inhomogeneous quadratic mapping
    d = len(x)
    r = np.zeros(d*(d+1)/2)
    k = 0
    for i in range(d):
        for j in range(i,d):
            if i==j: r[k] = x[i]*x[j]
            else: r[k] = np.sqrt(2)*x[i]*x[j]
            k += 1
    return r


kernel = "linear" #linear, quadratic, or Iquadratic (inhomogeneous)
loss = "quadratic" #quadratic or hinge

fname = sys.argv[1]
C = np.float(sys.argv[2])
if len(sys.argv) > 3:
    kernel = sys.argv[3]
if len(sys.argv) > 4:
    loss = sys.argv[4]

print "params", fname, C, kernel, loss

D = np.loadtxt(fname)
n,d = np.shape(D) #get input dimensions

X = D[:,0:d-1]
Y = D[:,d-1] #classes

if kernel == "linear":
    K = np.dot(X, X.T)+1 #have to add 1 for the w dim
    if loss == "quadratic":
        delta = np.diag(np.tile(1.0/(2*C), n))
        K = K+delta
elif kernel == "quadratic" or kernel == "Iquadratic":
    K = np.dot(X, X.T)
    if kernel == "Iquadratic":
        K = K+1 #add one to each element
    K = K*K #square each element
    if kernel == "quadratic":
        K = K+1 #add one for the last dimension, since we map to one
                #higher dimensionality for the bias
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

X = np.append(X, np.ones((n,1)), 1) #add extra column of ones
if kernel == "linear":
    ay = alpha*Y
    w = np.sum(X.T*ay, axis=1)
    print "hyperplane:", w
else: #quadratic
    if kernel == "Iquadratic":
        Z = np.array([phi2ih(x) for x in X]) #row-wise mapping to phi
    else:
	Z = np.array([phi2h(x) for x in X]) #row-wise mapping to phi
    ay = alpha*Y
    w = np.sum(Z.T*ay, axis=1)
    print "w:", w
    d = X.shape[1]
    print d
    if kernel == "Iquadratic":
        z = phi2ih(np.ones(d))
    else:
        z = phi2h(np.ones(d))
    print "h:", w*z


