#!/usr/bin/env python

import sys
import random
import numpy as np

kernel = "linear" #or quadratic
loss = "quadratic" #or hinge

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


w = np.tile(np.sqrt(C/d), d) 
w = np.zeros(d)
#start iterative stochastic gradient ascent
eps = 0.0001
t = 0
t0 = C
err = 1
while err > eps:
    w_prev = w.copy()
    fullXY = X.T*Y
    wXY = np.dot(w, fullXY)
    #print wXY
    idx = [i for (i,wxy) in enumerate(wXY) if wxy < 1]
    #print idx
    v = np.zeros(d)
    for i in idx:
        v += Y[i]*X[i,]

    if loss == "quadratic":
        S = np.zeros((d,d))
        for i in idx:
            xi = X[i,:].copy()
            xi.resize(d,1)
            S += np.dot(xi,xi.T)
        
        grad = w - 2*C*v + 2*C*np.dot(S,w)
        
        #compute best step size, eta
        #rat = 2*C*np.dot(np.dot(S,grad), grad)/np.dot(grad,grad)
        #eta = 1.0/(1.0 + rat)
        hess = np.eye(d) + 2*C*S
        HI = np.linalg.inv(hess)
        w = w - np.dot(HI,grad)
        #w = w - eta*grad
    else: #hinge
        grad = w - C*v
        eta = 1.0/(t0+t) 
        w = w - eta*grad

    t += 1
    err = np.linalg.norm(w-w_prev)
    print t, err

print "hyperplane:", w
