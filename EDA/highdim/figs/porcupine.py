#!/usr/bin/env python

import sys
from math import *
import numpy as np

def dfactorial(x): #double factorial, x should be odd
    if x <= 2: return x
    else: return x*dfactorial(x-2)

def volSphere(d):
    if d % 2 == 0:
        volS = pow(pi,d/2)/factorial(d/2)
    else:
        volS = pow(pi,(d-1)/2)*pow(2,(d+1)/2)/dfactorial(d)
    return volS

#def get_radius2(d):
    #if d == 2:
        #return 1.0/sqrt(2)
    #else:
        #volH = pow(2,d)
        #diff = volH-volSphere(d)
        #h = sqrt(diff*2*d/pow(2,d-1))
        #R=  1+h
        #return 1/R

def get_radius(d):
    volH = pow(2,d)
    diff = volH-volSphere(d)
    cd = sin(pi/pow(2,d-1))/sin(pi*(pow(2,d-1)-1)/pow(2,d))
    scd = sqrt(4-cd*cd)
    R = scd/2 + (diff + (pi - pow(2,d-2)*cd*scd))/(pow(2,d-1)*cd)
    return 1/R


####MAIN########

d = int(sys.argv[1]) #dimensionality
#r = float(sys.argv[2]) #inner rad
#R = float(sys.argv[3]) #outer rad

r = get_radius(d)
R = 1

N = pow(2,d)
incr = 360.0/N;
incr2 = incr/2.0;

print "%% d=%d r=%f R=%f" % (d, r, R)
print "\\psset{unit=1in,dotsep=0.5mm,PointName=none,PointSymbol=none}"
print "\\pstGeonode(0,0){O}"
print "\\pstGeonode(0,%0.3f){A}" % r
print "\\pstCircleOA[fillstyle=solid,fillcolor=lightgray]{O}{A}"

idx= 0
for a in np.arange(0,360,incr):
    ang = radians(a)
    point = (R*cos(ang), R*sin(ang))
    print "\\pstGeonode(%0.3f,%0.3f){R%d}" % (point[0], point[1], idx)
    idx += 1

for i in range(0,int(N)/2,1):
    a = i
    b = i+N/2
    print "\\pstLineAB[linestyle=dotted]{R%d}{R%d}" % (a,b)

idx= 0
for a in np.arange(incr2,360,incr):
    ang = radians(a)
    point = (r*cos(ang), r*sin(ang))
    print "\\pstGeonode(%0.3f,%0.3f){r%d}" % (point[0], point[1], idx)
    idx += 1

for n in range(0,int(N),1):
    a = n
    b = (n+1)%N
    print "\\pstLineAB{r%d}{R%d}" % (n,a)
    print "\\pstLineAB{r%d}{R%d}" % (n,b)


