#!/usr/bin/env python

import commands, sys
import multiprocessing as mp

def run_cmd(runnum):
    res=commands.getoutput(cmd)
    A = res.split()
    purity = 0.0
    dpurity = 0.0
    fval = 0.0
    for (j,w) in enumerate(A):
        if w == "Purity":
            purity = float(A[j+1])
            #print "purity", runnum, purity
        if w == "dPurity":
            dpurity = float(A[j+1])
            #print "dpurity", runnum, dpurity
        if w == "F":
            fval = float(A[j+1])
            print "F", runnum, fval


    return (purity, dpurity, fval, res)


#########
cmd = sys.argv[1]
numprocs = 8
n = 100
if len(sys.argv) > 2:
    n = int(sys.argv[2])


par = mp.Pool(numprocs)
runs = par.map(run_cmd, range(n))

maxpurity = 0.0
maxdpurity = 0.0
maxfval = 0.0
maxres = ""
i = 0
for (purity, dpurity, fval, res) in runs:
    if fval > maxfval: 
        maxfval = fval
        print i, maxfval
        maxres = res
    
    i += 1

print maxres
print "maxfval", maxfval
#print "maxpurity", maxpurity
#print "maxdpurity", maxdpurity

