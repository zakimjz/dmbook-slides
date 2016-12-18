#!/usr/bin/env python

f = open("iris-layout.txt", "rU")

E = []
P = []
V = {}
Vc = {}
Vd = {}

for l in f.readlines():
    a = l.strip().split()
    if a[0] == "e":
        E.append((a[1], a[2]))
    elif a[0] == "v":
        if a[4] not in Vc: Vc[a[4]] = []
        Vc[a[4]].append(int(a[1]))
        V[int(a[1])] = (a[2], a[3])
    elif a[0] == "d":
        for i in range(1,len(a)-1):
            P.append((a[i],a[i+1]))
        for i in range(1,len(a)):
             Vd[int(a[i])] = True

for n in sorted(V):
    print "\pnode(%s,%s){n%d}" % (V[n][0], V[n][1],n)


print "\psset{linecolor=lightgray}"
for e in E:
    print "\\ncline[]{n%s}{n%s}" % (e[0], e[1])

print "\psset{linecolor=black, linewidth=2pt}"
for e in P:
    print "\\ncline[]{n%s}{n%s}" % (e[0], e[1])


print "\psset{linecolor=black, linewidth=1pt}"
for cl in Vc:
    if cl == "Iris-setosa":
        print "\psset{dotstyle=Bsquare}"
    elif cl == "Iris-virginica":
        print "\psset{dotstyle=Btriangle}"
    else:
        print "\psset{dotstyle=Bo}"
    
    for n in sorted(Vc[cl]):
        if n in Vd:
            print "\dotnode[fillcolor=gray](%s,%s){n%d}"%(V[n][0], V[n][1],n)
        else:
            print "\dotnode[](%s,%s){n%d}"%(V[n][0], V[n][1],n)

