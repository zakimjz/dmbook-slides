#!/usr/bin/env python

import sys
from operator import itemgetter

f = open(sys.argv[1], "rU")
showclusters = int(sys.argv[2])
showcrossedges = int(sys.argv[3])
mcl = 1
if len(sys.argv) > 4:
    mcl = int(sys.argv[4])

if mcl:
    colary = ["gray","gray","lightgray"]
else:   
    colary = ["black","darkgray","lightgray"]


E = []
DE = []
P = []
V = {}
Vc = {}

for l in f.readlines():
    a = l.strip().split()
    if a[0] == "e":
        wt = float(a[3])
        dele = int(a[4])
        if wt >= 0.95: colidx = 0
        elif wt >= 0.9: colidx = 1
        else: colidx = 2
        E.append((a[1], a[2], colidx, dele))
    elif a[0] == "de":
        DE.append((a[1],a[2]))
    elif a[0] == "v":
        if a[5] not in Vc: Vc[a[5]] = []
        Vc[a[5]].append(int(a[1]))
        if len(a) > 6: att = int(a[6])
        else: att = 0
        V[int(a[1])] = (a[2], a[3], att)

#need to print edges first, otherwise, nodes will get written 
#over with lines, so print dummy nodes
for n in sorted(V):
    print "\pnode(%s,%s){n%d}" % (V[n][0], V[n][1],n)


print "\psset{linewidth=0.5pt,dotsep=2pt}"


for e in sorted(E, key=itemgetter(2), reverse=True):
    col = colary[e[2]]
    if (e[3] and showcrossedges):
        print "\\ncline[linestyle=dotted,linecolor=%s]{n%s}{n%s}"\
                    % (col, e[0], e[1])
    else:
        print "\\ncline[linecolor=%s]{n%s}{n%s}" % (col, e[0], e[1])

for e in DE:
    print "\\ncline[linecolor=black]{->}{n%s}{n%s}" % (e[0], e[1])

#print the fill circle nodes
for cl in Vc:
#    if cl == "Iris-setosa":
        #print "\psset{dotstyle=Bsquare}"
    #elif cl == "Iris-virginica":
        #print "\psset{dotstyle=Btriangle}"
    #else:
        #print "\psset{dotstyle=Bo}"

    if showclusters:
        if cl == "1": print "\psset{dotstyle=Bo}"
        elif cl == "2": print "\psset{dotstyle=Btriangle}"
        else: print "\psset{dotstyle=Bsquare}"
    else:
         print "\psset{dotstyle=Bo}"

    for n in sorted(Vc[cl]):
        if V[n][2]:
            print "\dotnode[fillcolor=gray](%s,%s){n%d}"%(V[n][0], V[n][1],n)
        else: 
            print "\dotnode[](%s,%s){n%d}"%(V[n][0], V[n][1],n)



