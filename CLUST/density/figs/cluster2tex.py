#!/usr/bin/env python

import sys
import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
from scipy import spatial
import multiprocessing as mp

#helper routines
clmap = {"Iris-versicolor": 0, "Iris-setosa": 1, 
         "Iris-virginica": 2}
def class2int(c):
	if dfname == "iris.txt":
		c = c.strip('"')
		return clmap[c]
	else:
		return int(c)

def density((xk,yk,myh)):
	val = 0.0
	vx = A[xk,yk]
	vy = B[xk,yk]
	eps = 3*myh
	Neps = kdtree.query_ball_point([vx, vy], eps)
	for i in Neps:
		x = map(float, X[i])
		a =	x[0] > 0 and "vx-%0.2f" % x[0] or "vx+%0.2f" % abs(x[0])
		b = x[1] > 0 and "vy-%0.2f" % x[1] or "vy+%0.2f" % abs(x[1])
		eqstr = "np.exp(%s*((%s)**2+(%s)**2))" % (gamma, a, b)
		#print eqstr
		val += eval(eqstr)
	val *= (1.0/(n*myh*np.sqrt(2*np.pi)))
	print xk, yk, vx, vy, val
	return val



#######################

xi = 0.0
dname = sys.argv[1]
if len(sys.argv) > 2:
    xi = float(sys.argv[2]) #use to filter the density vals

if dname == "iris":
    dfname = "iris.txt" #data file
    #cfname = "iris-2d-h0.2.out" #cluster output file
    #ofname = "iris-clustered.out"
    ofnameps = "iris-pspoints.tex"
    classattr = 4
    xyscale = 1.0
    zscale = 10.0
    hvals2d = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5]
    hvals1d = [0.25, 0.5, 1.0, 2.0]
    rangestep = (4,8,0.04)
    A = np.arange(3.5, 8.6, 0.075)
    B = np.arange(1.0, 5.1, 0.075)
    hval2d = [0.2]
elif dname == "t7-4k":
    dfname = "t7-4k.txt" #data file
    #cfname = "t7-4k-hX.out" #cluster output file
    #ofname = "t7-4k-clustered.out"
    ofnameps = "t7-4k-pspoints.tex"
    classattr = 2
    xyscale = 100.0
    zscale = 20000.0
    hvals2d = [5, 10, 15, 20, 25, 30, 35, 40]
    hvals1d = [15, 20, 25, 30]
    rangestep = (20,600,10)
    A = np.arange(0.0, 705.0, 7.5)
    B = np.arange(0.0, 505.0, 7.5)
    hval2d = [10]

df = open(dfname, "rU")
#cf = open(cfname, "rU")
#of = open(ofname, "w")
ofps = open(ofnameps, "w")

X = []
Y = []
n = 0
for l in df.readlines():
    a = l.strip().split(",")
    pt = a[:2] #get only the first 2 dims
    pta = ["%0.4f"%float(x) for x in pt]
    X.append(pta)
    Y.append(class2int(a[-1]))
    print >>ofps, "\\psPoint(%s, %s, %f){p%d}" %\
            (float(pta[0])/xyscale, float(pta[1])/xyscale, xi, n)
    print >>ofps, "\\psdot(p%d)" % (n)
    n += 1
ofps.close()
print n

#2d Gaussian Density function
gamma = "(-0.5/\myh^2)"
equation = []
for i in range(n):
    x = map(float, X[i])
    a = x[0] > 0 and "x-%0.2f" % x[0] or "x+%0.2f" % abs(x[0])
    b = x[1] > 0 and "y-%0.2f" % x[1] or "y+%0.2f" % abs(x[1])
    eqstr = "e^(%s*((%s)^2 + (%s)^2))" % (gamma, a, b)
    equation.append(eqstr)

ff1 = "(%0.6f/\myh^2)" % (1.0/(n*2*np.pi))
part2 = "+\n".join(equation)
finalstr =  "%s*\n(%s)" % (ff1, part2)
print finalstr


#shortcut the influence
D = np.loadtxt(dfname,delimiter=",", 
               usecols = (0,1,classattr),
               converters={classattr:class2int})

kdtree = spatial.KDTree(D[:,0:2])

#fig = plt.figure()
#ax = fig.gca(projection='3d')
gamma = "(-0.5/myh**2)"
nA = len(A)
nB = len(B)
totN = nA*nB
A, B = np.meshgrid(A, B)
A = A.T
B = B.T
xn, yn = A.shape
C = A*0.0
E = A*0.0
tups = []
eps = 0.0
Neps = []
numprocs = 10
par = mp.Pool(numprocs)
for myh in hval2d:
    fsname = dname+"-h"+str(myh)+".surface.obj"
    fs = open(fsname, "w")
    maxval = 0.0
    for xk in range(xn):
        for yk in range(yn):
            tups.append((xk, yk, myh))
    vals = par.map(density, tups)

    cc = 0
    for xk in range(xn):
        for yk in range(yn):
            val = vals[cc]
            if val > maxval: maxval = val
            C[xk,yk] = val*zscale
            #E[xk,yk] = xi
            cc += 1
            print xk, yk, C[xk,yk]
            if C[xk,yk] <= xi:
                C[xk,yk] = xi
            print >>fs, "v %0.5f %0.5f %0.5f" %\
                (A[xk,yk]/xyscale, B[xk,yk]/xyscale, C[xk,yk])
            #ax.scatter(A[xk,yk],B[xk,yk], val, c='r', marker='.', s=10)

    print "myh", myh, maxval

    print A
    print B
    print C
    
    f1 = []
    f2 = []
    f3 = []
    f4 = []
    for i in range(2,totN-nB,nB):
        for j in range(1, nB):
            ii = i+j-1
            jj = ii-1
            f1.append(ii)
            f2.append(jj)
            f3.append(jj+nB)
            f4.append(jj+nB+1)
    for i in range(2,totN-nB,nB):
        for j in range(1, nB):
            ii = i+j-1
            jj = ii-1
            f4.append(ii)
            f3.append(jj)
            f2.append(jj+nB)
            f1.append(jj+nB+1)

    for a in range(len(f1)):
        print >>fs, "f", f1[a], f2[a], f3[a], f4[a] 
    fs.close()
    #ax.plot_surface(A, B, E, rstride=1, cstride=1, linewidth=0, alpha=0.85)
    #ax.plot_surface(A, B, C, rstride=1, cstride=1, linewidth=0,
    #		cmap=cm.jet, alpha=0.75)
    #ax.plot_wireframe(A, B, C, rstride=1, cstride=1, linewidth=1)
    #cset = ax.contour(A, B, C, zdir='z', offset=-0.2)
    #plt.savefig('surface.eps')
    #plt.show()


if dname == "iris":
    #1d Gaussian function: based on X[0]
    equation = []
    gamma = "(-0.5/\myh^2)"
    for i in range(n):
        x = map(float, X[i])
        a = x[0] > 0 and "x-%0.2f" % x[0] or "x+%0.2f" % abs(x[0])
        eqstr = "Euler^(%s*((%s)^2))" % (gamma, a)
        equation.append(eqstr)
    ff1 = "(%0.6f/\myh)" % (1.0/(n*np.sqrt(2*np.pi)))
    print "%s*\n(%s)" % (ff1, "+\n".join(equation))


    #gamma = "(-0.5/myh**2)"
    #for myh in hvals2d:
    #    maxval = 0
    #    for vx in np.arange(rangestep[0], rangestep[1], rangestep[2]):
    #        val = 0
    #        for i in range(n):
    #            x = map(float, X[i])
    #            a = x[0] > 0 and "vx-%0.2f" % x[0] or "vx+%0.2f" % abs(x[0])
    #            eqstr = "np.exp(%s*((%s)**2))" % (gamma, a)
    #            #print eqstr
    #            val += eval(eqstr)
    #        val *= (1.0/(n*myh*np.sqrt(2*np.pi)))
    #        if val > maxval: maxval = val
    #    print myh, maxval


    #1d Discrete function: based on X[0]
    for h in hvals1d:
        f = open(dname+"1d-h%0.2f.txt"%h, "w")
        maxval = 0
        for vx in np.arange(rangestep[0], rangestep[1], rangestep[2]):
            val = 0
            for i in range(n):
                x = map(float, X[i])
                if abs((vx-x[0])/h) <= 0.5:
                    val += 1
            val *= 1.0/(n*h)
            if val > maxval: maxval = val
            print >>f, vx, val
        f.close()

        print h, maxval


    #C = {}
    #mu = {}
    #N = {}
    #for l in cf.readlines():
    #    a = l.strip().translate(None,"[],").split() #delete all [],
    #    if a[0] == "cluster":
    #        cid = a[1]
    #    elif a[0] == "mean":
    #        mu[cid] = " ".join(a[1:])
    #    elif a[0] == "pts":
    #        C[cid] = [int(x) for x in a[1:]]
    #    elif a[0] == "size":
    #        N[cid] = a[1]

    #for i in sorted(C):
    #    #print "cluster", i
    #    #print N[i], "--", mu[i]
    #    #print len(C[i]), C[i]
    #    TC = {}
    #    for j in C[i]:
    #        if Y[j] not in TC: TC[Y[j]] = []
    #        TC[Y[j]].append(j)
    #    for id in TC:
    #        print >>of, "%", i, id, len(TC[id])
    #        for j in TC[id]:
    #            x = map(float, X[j])
    #            print >>of, str(x).translate(None,"[]")+", 0"


#part12 = "2.7183 x dup mul y dup mul add %0.2f neg mul exp %f mul"%\
        #(abs(gamma), ff1)
#part22 = " add %\n".join(equation2[1:])
#part22 = " %\n"+equation2[0]+"% \n"+part22+" add "
#print "%s %s mul" % (part12, part22)

#h=0.3
#gamma = -0.5/h**2
#sval = 0
#for vx in np.arange(-4,4,0.5):
    #for vy in np.arange(-2,2,0.5):
        #print vx, vy
        #for i in range(150):
            #x = map(float, X[i])
            #expstr = "-2*(%0.2f*vx+%0.2f*vy)+%0.2f" %\
            #(x[0], x[1], x[0]**2+x[1]**2)
            #eqstr = "np.exp(%f*(%s))" % (gamma, expstr)
            ##print eqstr
            #val = eval(eqstr)
            #sval += val
        #print sval
