#!/usr/bin/env python

import numpy as np
import numpy.linalg as la


X = np.loadtxt('iris-slwpl.dat')
#print X
m = np.mean(X, axis=0)
#print m

Z = X-m

#save the centered data
np.savetxt('iris-slwpl-c.dat', Z, fmt="%0.4f")

a = np.min(Z, axis=0)
b = np.max(Z, axis=0)
print a
print b
print b-a

print "cov matrix"
nr, nc = np.shape(Z)
S = np.cov(np.transpose(Z))*(nr-1)/nr
print S

w, v = la.eigh(S)
print "eigvals, eigvecs"
print w
print v

#verify
SS = np.zeros((nc,nc))
for i in range(nc):
    vv = np.asmatrix(v[:,i])
    SS = SS + w[i]*np.dot(vv.T, vv)
print "SS", SS

#find sum of squared magnitudes
SS = 0
for i in range(nr):
    SS += la.norm(Z[i,:])**2
print "squared sum", SS

#find projected points
U = np.asmatrix([v[:,2],v[:,1],v[:,0]]).T
print "U", U
P = np.dot(U.T, Z.T).T
#print P
#save the projected data
np.savetxt('iris-slwpl-p.dat', P, fmt="%0.4f")
f = open("iris-slwpl-p-dots.tex", "w")
scale=np.array([1, 0.4, 0.2])
for i in range(nr):
    pt1 = np.squeeze(P[i,:].A)
    pt1 = pt1*scale
    pt1 = ["%0.3f" % x for x in pt1]
    pt1 = ", ".join(pt1)
    print >>f, "\\psPoint(%s){p%d}" % (pt1, i)
    print >>f, "\\psdot(p%d)" % i
f.close()

a = np.min(P, axis=0)
b = np.max(P, axis=0)
print a
print b
print b-a
c = b-a
C = np.diag(np.squeeze(c.A))
U*C
A = np.diag(np.squeeze(a.A))
B = np.diag(np.squeeze(b.A))
U*A
U*B

A2 = A
A2[2,2]=0
np.sum(U*A2, axis=1) #origin of PC1+PC2

A3 = A
A3[0,0]=0
np.sum(U*A3, axis=1)

print "scaled eigvecs, ui*lamdai"
print w[0]*v
print w[1]*v
print w[2]*v

print "check that eigvecs are orthogonal"
D = np.dot(np.transpose(v), v)
print D

#project points onto the principal component
u1 = np.asmatrix(v[:,2])
print u1
Q = np.asarray(np.dot(u1.T,u1))
print "projection matrix"
print Q

print "mapped points"
Y = np.dot(Z,Q)
#save the centered data
np.savetxt('iris-1dproj.dat', Y, fmt="%0.4f")

f1d = open('iris-1dmap.tex', "w")

for i in range(nr):
    pt1 = ["%0.3f" % x for x in Z[i,:]]
    pt1 = ", ".join(pt1)
    pt2 = ["%0.3f" % x for x in Y[i,:]]
    pt2 = ", ".join(pt2)
    print >>f1d, "\\pstThreeDLine[linecolor=gray]"+\
    "("+pt1+")"+"("+pt2+")"


#eigen sorts eig vals in increasing order
#project points onto two principal components
U = np.asmatrix([v[:,2],v[:,1]]).T
print "U matrix"
print U
Q2 = np.asarray(np.dot(U,U.T))
print "projection matrix Q2"
print Q2

print "mapped points"
Y2 = np.dot(Z,Q2)
#save the projected data
np.savetxt('iris-2dproj.dat', Y2, fmt="%0.4f")
A2 = np.dot(U.T,Z.T).T
#save reduced dim data
np.savetxt('iris-2dproj-dimr.dat', A2, fmt="%0.4f")
print "A mat", min(A2[:,0]), max(A2[:,0])
print "A mat", min(A2[:,1]), max(A2[:,1])


#project points onto non optimal axes
Uno = np.asmatrix([v[:,1],v[:,0]]).T
print "Uno matrix"
print Uno
Q2no = np.asarray(np.dot(Uno,Uno.T))
print "projection matrix Q2no"
print Q2no

print "mapped points NO"
Y2no = np.dot(Z,Q2no)



f2d = open('iris-2dmap.tex', "w")
fL = open('iris-2dL.dat', "w") #halfplane points
fR = open('iris-2dR.dat', "w")
f2dno = open('iris-2dmapnonopt.tex', "w")
fLno = open('iris-2dLnonopt.dat', "w") #halfplane points
fRno = open('iris-2dRnonopt.dat', "w")

for i in range(nr):
    pt = ["%0.3f" % x for x in Z[i,:]]
    pt1 = ", ".join(pt)
    norm= np.array([-0.79, 0.79, 0.42])
    if np.dot(norm,Z[i,:]) >= 0:
        hp = "L"
    else:
        hp = "R"

    if hp == "L":
        print >>fL, " ".join(pt)
    else:
        print >>fR, " ".join(pt)

    nonoptnorm = np.array([-1.43, 0.32, -3.36])
    if np.dot(nonoptnorm,Z[i,:]) >= 0:
        nonopthp = "L"
    else:
        nonopthp = "R"
    if nonopthp == "L":
        print >>fLno, " ".join(pt)
    else:
        print >>fRno, " ".join(pt)


for i in range(nr):
    pt1 = ["%0.3f" % x for x in Z[i,:]]
    pt1 = ", ".join(pt1)
    norm= np.array([-0.79, 0.79, 0.42])
    if np.dot(norm,Z[i,:]) >= 0:
        lcoll="black"
    else:
        lcoll = "gray"
    pt2 = ["%0.3f" % x for x in Y2[i,:]]
    pt2 = ", ".join(pt2)
    print >>f2d, "\\pstThreeDLine[linecolor="+lcoll+"]"+\
    "("+pt1+")"+"("+pt2+")"
f2d.close()

for i in range(nr):
    pt1 = ["%0.3f" % x for x in Z[i,:]]
    pt1 = ", ".join(pt1)
    norm= np.array([-1.43, 0.32, -3.36])
    if np.dot(norm,Z[i,:]) >= 0:
        lcoll="gray"
    else:
        lcoll = "black"
    pt2 = ["%0.3f" % x for x in Y2no[i,:]]
    pt2 = ", ".join(pt2)
    print >>f2dno, "\\pstThreeDLine[linecolor="+lcoll+"]"+\
    "("+pt1+")"+"("+pt2+")"
f2dno.close()



#f2d = open('iris-2dPCA.tex', "w")
#for i in range(nr):
    #pt1 = ["%0.3f" % (x,) for x in Z[i,:]]
    #pt1 = ", ".join(pt1)
    #ptid = "p"+str(i)
    #norm= np.array([-0.79, 0.79, 0.42])
    #if np.dot(norm,Z[i,:]) >= 0:
        #coll = "white"
    #else:
        #coll = "gray"
    #print >>f2d, "\\psPoint("+pt1+"){"+ptid+"}"
    #print >>f2d, "\psdot[fillcolor="+coll+"]("+ptid+")"
#f2d.close()

#SVD
L, S, R = la.svd(Z, 0)
print "SVD"
print S
print L
print R.T
