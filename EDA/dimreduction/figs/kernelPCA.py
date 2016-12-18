#!/usr/bin/env python

import numpy as np
import numpy.linalg as la


def Phi(Z):
    pmat = np.zeros((n,6))
    for i in range(n):
        pmat[i,:] = np.matrix([1, np.sqrt(2)*Z[i,0],
                               np.sqrt(2)*Z[i,1],
                               np.sqrt(2)*Z[i,0]*Z[i,1],
                               Z[i,0]**2, Z[i,1]**2])

    return pmat

def PhiH(Z):
    pmat = np.zeros((n,3))
    for i in range(n):
        pmat[i,:] = np.matrix([np.sqrt(2)*Z[i,0]*Z[i,1],
                               Z[i,0]**2, Z[i,1]**2])

    return pmat


def projectOnPC(a,ipc):
    if Ktype == "linear":
        print a.T*Z[:,0]
        print a.T*Z[:,1]
    elif Ktype == "quadratic":
        print np.sum(a)
        print 2*a.T*Z[:,0]
        print 2*a.T*Z[:,1]
        print 2*a.T*(Z[:,0].A*Z[:,1].A) #use array/elementwise
                    #multiplication for the two columns of Z
        print 2*a.T*(Z[:,0].A*Z[:,0].A)
        print 2*a.T*(Z[:,1].A*Z[:,1].A)
    elif Ktype == "quadraticH":
        print 2*a.T*(Z[:,0].A*Z[:,1].A) #use array/elementwise
                    #multiplication for the two columns of Z
        print 2*a.T*(Z[:,0].A*Z[:,0].A)
        print 2*a.T*(Z[:,1].A*Z[:,1].A)
    elif Ktype == "gaussian":
        xixi = np.mat(np.exp(-gamma*np.diag(Z*Z.T))).T
        print np.shape(xixi), np.shape(a)
        c = a.A*xixi.A
#        print "\\mye^(-%0.3f*(x^2+y^2)) * ( 0 "%(gamma)
        #for i in range(n):
            #if c[i] >= 0.00001:
                #print "+ %0.5f*\\mye^(%0.3f*(%0.3f*x+%0.3f*y)) "\
                        #% (c[i], 2*gamma, Z[i,0],Z[i,1])
        #print ")"
        


#Ktype = "linear"
#Ktype = "quadratic"
Ktype = "quadraticH" #homogeneous
#Ktype = "gaussian"
sigma = 0.01
gamma = 1/(2*sigma**2)

Z = np.matrix(np.loadtxt('iris-2d-nonlinear.dat'))
n, d = np.shape(Z)

#create kernel matrix
K = np.mat(np.zeros((n,n)))

for i in range(n):
    for j in range(i,n):
        if Ktype == "quadratic":
            kern = (1 + Z[i,:]*Z[j,:].T)**2 #quadratic kernel
        elif Ktype == "quadraticH":
            kern = (Z[i,:]*Z[j,:].T)**2 #homogeneous, quadratic
        elif Ktype == "linear":
            kern = Z[i,:]*Z[j,:].T
        elif Ktype == "gaussian":
            kern = np.exp(-gamma*la.norm(Z[i,:]-Z[j,:])**2)

        K[i,j] = kern
        K[j,i] = kern


#center in kernel space
I = np.mat(np.eye(n))
On = np.mat(np.ones((n,n))/n)


Kc = (I-On)*K*(I-On)

kw, kv = la.eigh(Kc)
print "eigvals, eigvecs"
print kw[n-10:]
print kw[n-10:]/n

#scale the eigenvector by the sqrt of eigenval
a1 = np.asmatrix(kv[:,n-1]*np.sqrt(1/kw[n-1]))
#print a1
#print la.norm(a1)**2

#project each point onto kernel PC1
A1 = Kc*a1
#print A1

#project each point onto kernel PC2
a2 = np.asmatrix(kv[:,n-2]*np.sqrt(1/kw[n-2]))
A2 = Kc*a2

#project each point onto kernel PC3
a3 = np.asmatrix(kv[:,n-3]*np.sqrt(1/kw[n-3]))
A3 = Kc*a3

A = np.hstack([A1,A2,A3])

print "lambdas", kw[n-1]/n, kw[n-2]/n, kw[n-3]/n
print kw[n-10:n]/n
print np.sum(kw/n)
print np.sum(kw[n-2:]/n)/np.sum(kw/n)

np.savetxt('iris-2d-kPC-'+Ktype+'.dat', A[:,0:2], fmt="%0.4f")
print "projected coordinates"
print np.min(A[:,0]), np.max(A[:,0])
print np.min(A[:,1]), np.max(A[:,1])
print np.min(A[:,2]), np.max(A[:,2])
print np.mean(A[:,0]), np.mean(A[:,1]), np.mean(A[:,2])

print "cov", np.cov(A.T)*(n-1.0)/n

#recover the principal components
if Ktype == "linear":
    PZ = Z
elif Ktype == "quadratic":
    #inhomogeneous quadratic in 2D
    PZ = Phi(Z)
elif Ktype == "quadraticH":
    PZ = PhiH(Z)

if Ktype != "gaussian":
    u1 = PZ.T*a1
    print u1
    print la.norm(u1)
    u2 = PZ.T*a2
    print u2

#lines of constant orthogonal projection

print "project on a1"
projectOnPC(a1,1)
print "project on a2"
projectOnPC(a2,2)
print "project on a3"
projectOnPC(a3,3)

#verify the results
if Ktype != "gaussian":
    m = np.mean(PZ, axis=0)
    PZc = PZ-m
    S = np.cov(np.transpose(PZc))*(n-1.0)/n    
    print S
    w, v = la.eigh(S)
    print "eigvals, eigvecs"
    print w
    print v
