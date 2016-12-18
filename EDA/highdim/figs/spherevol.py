#!/usr/bin/env python

from math import *

d=50
for i in range(1,d+1):
    vd = pow(sqrt(pi),i)/gamma((i+2)/2.0)
    print i, vd

