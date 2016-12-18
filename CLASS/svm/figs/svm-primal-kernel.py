#!/usr/bin/env python

import sys
import random
import numpy as np

def phi2h(x): #inhomogeneous quadratic mapping
    r = np.zeros(d*(d+1)/2)
    k = 0
    for i in range(d):
        for j in range(i,d):
            if i==j: r[k] = x[i]*x[j]
            else: r[k] = np.sqrt(2)*x[i]*x[j]
            k += 1
    return r


kernel = "linear" #default
loss = "quadratic" #default, versus hinge

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
elif kernel == "quadratic":
    K = np.dot(X, X.T)
    K = K*K #square each element

b = np.zeros(n) #initialize b

#start iterative stochastic gradient ascent
eps = 0.0001
t = 0
t0 = C
ridge = 1
err = 1

while err > eps:
    b_prev = b.copy()

    fullKY = K.T*Y
    bKY = np.dot(b, fullKY)
    #print bKY
    idx = [i for (i,bky) in enumerate(bKY) if bky < 1]
    #print idx
    v = np.zeros(n)
    for i in idx:
        v += Y[i]*K[i,:]

    if loss == "quadratic":
        S = np.zeros((n,n))
        for i in idx:
            ki = K[i,:].copy()
            ki.resize(n,1)
            S += np.dot(ki, ki.T)

        K2CS = K+2*C*S
        grad = np.dot(K2CS,b) - 2*C*v 
        #eta = np.dot(grad,grad)/np.dot(grad, np.dot(K2CS,grad))

        hess = K2CS #hess is the same as K2CS
        hess = K2CS + np.diag(np.tile(ridge,n))
        #b = b - eta*grad
        eta = 1
    else: #loss == hinge
        grad = np.dot(K,b) - C*v
        hess = K+ np.diag(np.tile(ridge,n))
        eta = 1.0/(t0+t)

    HI = np.linalg.inv(hess) #inverse of hessian
    b = b - eta*np.dot(HI, grad)

    t += 1
    err = np.linalg.norm(b-b_prev)
    print t, err

print "coefficients:", b

if kernel == "linear":
    w = np.sum(X.T*b, axis=1)
    print "hyperplane:", w
else: #quadratic inhomogeneous
    Z = np.array([phi2h(x) for x in X]) #row-wise mapping to phi
    w = np.sum(Z.T*b, axis=1)
    print "w:", w
    z = phi2h(np.ones(d))
    print "h:", w*z


