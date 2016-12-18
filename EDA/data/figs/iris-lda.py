#!/usr/bin/env python

Vc = {}

f = open("iris-slw-clbl.txt", "rU")
for l in f.readlines():
    a = l.strip().split()
    if a[0] not in Vc: Vc[a[0]] = []
    Vc[a[0]].append((float(a[1]),float(a[2])))

#print "\pstGeonode[PointSymbol=none](1.0738,-1.3766){na}"
#print "\pstGeonode[PointSymbol=none](-1.0738,1.3766){nb}"

n=1;
for cl in Vc:
    if cl == "Iris-setosa": 
        ps = "Bsquare"
        pps = "square"
    elif cl == "Iris-virginica": 
        ps = "Btriangle"
        pps = "triangle"
    else: 
        ps = "Bo"
        pps = "o"
    
    for (a,b) in Vc[cl]:
        print "\pstGeonode[PointSymbol=%s,fillcolor=lightgray](%0.2f,%0.2f){n%d}" % (ps,a,b,n)
        print "\pstProjection[PointSymbol=%s,dotscale=1.25]{na}{nb}{n%d}[p%d]" % (pps,n,n)
        print "\pstLineAB[linestyle=dotted,linecolor=gray]{n%d}{p%d}" % (n,n)
        n = n+1
