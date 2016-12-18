#!/usr/bin/env python
import sys
RL = {}
for fn in sys.argv[1:]:
    f = open(fn,"rU")
    for l in f:
        a = l.strip().split(';')
        vals = a[0].strip().split()
        iset = a[1].strip()
        if iset not in RL: RL[iset] = []
        lorig = float(vals[8])
        lswap = float(vals[9])
        RL[iset].append((lorig-lswap)/lorig)

ARL = []
for iset in RL:
    avgrl = sum(RL[iset])/len(RL[iset])
    ARL.append((avgrl, iset))

for (r,i) in sorted(ARL):
    print r, i 

#iset = "11 42"
#print iset, " -- ", avgrl
#PMF = {}
#for rl in RL[iset]:
#    if rl not in PMF: PMF[rl] = 0.0
#    PMF[rl] += 1
#
#s = 0
#for rl in sorted(PMF):
#    print rl, PMF[rl]/len(RL[iset])
#    s += PMF[rl]/len(RL[iset])
#print s
