#!/usr/bin/env python

import sys
import networkx as nx


def read_graphs(fname):
    #read graphs
    f = open(fname, "rU")
    D = {}
    I = {}
    E = {}
    for l in f.readlines():
        a = l.strip().split()
        if a[0] == 't': #graph id 
            gid = int(a[2])
            G = nx.Graph(gid=gid)
            D[gid] = G
        elif a[0] == 'v': #vertices
            vid = int(a[1])
            vlbl = a[2]
            G.add_node(vid, L=vlbl)
            if vlbl not in I: I[vlbl] = {}
            I[vlbl][gid] = True
        elif a[0] == 'e': #edges
            vi = int(a[1])
            vj = int(a[2])
            elbl = 'e'
            if len(a) > 3: 
                elbl = a[3]
            li = G.node[vi]['L']
            lj = G.node[vj]['L']
            if li <= lj:
                edge = (0, 1, li, lj, elbl)
            else:
                edge = (0, 1, lj, li, elbl)
            G.add_edge(vi, vj, L=elbl)
            prefix = edge[2]
            if prefix not in E: E[prefix] = {}
            if edge not in E[prefix]: E[prefix][edge] = {}
            E[prefix][edge][gid] = True
    
    #for G in D:
    #    print G.graph['gid']
    #   print G.nodes(), G.nodes(data=True)
    #    print G.edges()

    return D, I, E

def edgecompare((vi, vj), (vx, vy)):
    fij = True
    fxy = True
    if vi > vj: fij = False
    if vx > vy: fxy = False

    if fij and fxy: #both are forward
        ff = cmp(vj, vy)
        if ff == 0: ff = cmp(vx,vi)
        return ff
    elif not fij and not fxy: #both are backward
        bb = cmp(vi, vx)
        if bb == 0: bb = cmp(vj,vy)
        return bb
    elif fij and not fxy: #forward, backward
        fb = 1
        if vj <= vx: fb = -1
        return fb
    else: #backward, forward
        bf = 1
        if vi < vy: bf = -1
        return bf

def dfscompare((vi, vj, li, lj, eij), (vx, vy, lx, ly, exy)):
    #print vi, vj, li, lj, vx, vy, lx, ly
    ecmp = edgecompare((vi,vj), (vx,vy))
    if ecmp == 0: #compare labels
        return cmp((li, lj, eij), (lx, ly, exy))
    else:
        return ecmp


#return the nodes on the rightmost path
def get_rmp(C): 
    rp = [0] #nodes on the rightmost path
    rc = 0 #rightmost child (last element of rp)
    if len(C) > 0:
        for t in C:
            (vi, vj, li, lj, eij) = t
            if vj > vi:
                v = rp.pop()
                while(v != vi):
                    v = rp.pop()
                rp.append(vi)
                rp.append(vj)
        rc = rp[-1]

    return (rp, rc)


def subgraph_isomorphisms(C, G):
    (vi, vj, li, lj, eij) = C[0] #first tuple
    L = nx.get_node_attributes(G,'L')
    #find nodes with label li
    M = [[x] for x in L if L[x]==li] 
    #print "lset", gid, L, M

    #process remaining tuples
    for (vi, vj, li, lj, eij) in C:
        #print "now do", vi, vj, li, lj
        new_M = []
        for f in M:
            fv = f[vi]
            if vj > vi:#forward extensions
                valid = [x for x in G.neighbors(fv)
                            if x not in f and
                               G.node[x]['L'] == lj and
                               G.edge[fv][x]['L'] == eij]
                for x in valid:
                    new_f = list(f)
                    new_f.append(x)
                    new_M.append(new_f)
            else: #backward extension
                fx = f[vj]
                if fx in G.neighbors(fv):
                    new_M.append(f)
        #update the isomporphisms
        M = new_M

    return M


def get_isomorphisms(C, DB):
    #print "xx", C
    M = {}
    for gid in DB:
        G = DB[gid]
        if len(C) > 0: #candidate is not empty
            #get all isomorphisms
            ISO = subgraph_isomorphisms(C, G)
            if len(ISO) > 0:
                M[gid] = ISO
    return M


# return set of rightmost path extensions
def frmpe(C, DB, leastext=False):
    #print "rmpe", C
    rmp, rmc = get_rmp(C)
    #print rmp, rmc
    E = {}
    for gid in DB:
        G = DB[gid]
        M = []
        if len(C) > 0: #candidate is not empty
            #get all isomorphisms
            M = subgraph_isomorphisms(C, G)
        
            #get all extensions
            for f in M:
                iso = [(v,f[v]) for v in rmp]
                iso.reverse() #process from rmc to the root
                #print "M", f, iso
                
                #backward extension
                vf = {}
                for v in rmp: 
                    vf[f[v]] = v #inverse iso from rmp
                (v, fv) = iso[0] #rmc mapping
                for x in G.neighbors(fv):
                    if x in vf:
                        e = (v,vf[x],G.node[fv]['L'], G.node[x]['L'],\
                             G.edge[fv][x]['L'])
                        be = (vf[x], v, G.node[x]['L'], G.node[fv]['L'],\
                              G.edge[x][fv]['L'])
                        if e not in C and be not in C:
                            if e not in E: E[e] = {}
                            E[e][gid] = True
                            #print "added back", fv, x, e

                #forward extensions
                for v, fv in iso:
                    #print v, fv
                    for x in G.neighbors(fv):
                        if x not in f:
                            e = (v, rmc+1, G.node[fv]['L'], G.node[x]['L'],\
                                 G.edge[fv][x]['L'])
                            if e not in E: E[e] = {}
                            E[e][gid] = True 
                            #print "added for", fv, x, e
        else: #candidate is empty
            #add all edges in G into E
            for (u,v) in G.edges():
                e = (0,1,G.node[u]['L'], G.node[v]['L'],
                    G.edge[u][v]['L'])
                if e not in E: E[e] = {}
                E[e][gid] = True
                #reverse the nodes
                if G.node[v]['L'] != G.node[u]['L']:
                    e = (0,1,G.node[v]['L'], G.node[u]['L'],
                        G.edge[v][u]['L'])
                    if e not in E: E[e] = {} 
                    E[e][gid] = True
               
    Eset = []
    for e in sorted(E, cmp=dfscompare):
        s = len(E[e])
        if not leastext and s < minsup: continue
        Eset.append((e, s))

    return Eset

# return set of rightmost path extensions
def rc_rmpe(C, ISO, DB, leastext=False):
    #print "rmpe", C
    rmp, rmc = get_rmp(C)
    #print rmp, rmc
    E = {}
    for gid in ISO:
        G = DB[gid]
        #get all isomorphisms
        M = ISO[gid]
        
        #get all extensions
        for f in M:
            iso = [(v,f[v]) for v in rmp]
            iso.reverse() #process from rmc to the root
            #print "M", f, iso
            
            #extensions from rc
            vf = {}
            for v in rmp: 
                vf[f[v]] = v #inverse iso from rmp
            (v, fv) = iso[0] #rmc mapping
            for x in G.neighbors(fv):
                #backward
                if x in vf:
                    e = (v,vf[x],G.node[fv]['L'], G.node[x]['L'],\
                         G.edge[fv][x]['L'])
                    be = (vf[x], v, G.node[x]['L'], G.node[fv]['L'],\
                          G.edge[x][fv]['L'])
                    if e not in C and be not in C:
                        if e not in E: E[e] = {}
                        E[e][gid] = True
                        #print "added back", fv, x, e

                #forward extensions
                if x not in f:
                    e = (v, rmc+1, G.node[fv]['L'], G.node[x]['L'],\
                         G.edge[fv][x]['L'])
                    if e not in E: E[e] = {}
                    E[e][gid] = True 
                    #print "added for", fv, x, e
               
    Eset = []
    for e in sorted(E, cmp=dfscompare):
        s = len(E[e])
        if not leastext and s < minsup: continue
        Eset.append((e, s))

    return Eset


#convert dfscode into graph
def make_graph(C):
    G = nx.Graph(gid=0)
    for (vi, vj, li, lj, eij) in C:
        G.add_node(vi, L=li)
        G.add_node(vj, L=lj)
        G.add_edge(vi, vj, L=eij)
    return G
        

def iscanonical(C):
    G = make_graph(C)
    #print "\niscan", C
    #print G.graph, G.edges(), G.nodes()
    DB = {0: G}
    nC = []
    for i in range(len(C)):
        #print "cann", nC
        E = frmpe(nC, DB, leastext=True)
        #print E
        e,s = E[0] #look at the least extension
        #print "compare", e, i, C[i]
        if dfscompare(e,C[i]) == -1:
            return False
        nC.append(e)

    return True

def gspan(C, M, DB, EXT):
    gspan.cnt += 1 #static graph counter
    rmp, rmc = get_rmp(C)
    print gspan.cnt, C
    E = rc_rmpe(C,M,DB) #frequent rightmost path extensions & their supports
    #print "rcext", E
    nEXT = []
    for (e,s) in EXT:
        ne = (e[0], rmc+1, e[2], e[3], e[4])
        nEXT.append((ne,s))

    #print "ext", nEXT
    E.extend(nEXT)
    #print "all", E
    for (i, ext) in enumerate(E):
        (e,s) = ext
        nC = list(C)
        nC.append(e)
        nM = get_isomorphisms(nC, DB)
        if iscanonical(nC) and len(nM) >= minsup:
            nEXT = E[i:]
            gspan(nC, nM, DB, nEXT)

######## main #########
fname = sys.argv[1] #graph filename
minsup = int(sys.argv[2]) #absolute min sup 

gspan.cnt = 0
DB, V, E = read_graphs(fname)

for l in V:
    sup = len(V[l])
    if sup >= minsup:
        print l, sup

EXT = {}
for a in sorted(E):
    for e in sorted(E[a], cmp=dfscompare):
        sup = len(E[a][e])
        if sup >= minsup:
            if a not in EXT: EXT[a] = []
            EXT[a].append((e, sup))

for a in EXT:
    for i in range(len(EXT[a])):
        C,s = EXT[a][i]
        extns = EXT[a][i:]
        M = get_isomorphisms([C], DB)
        gspan([C], M, DB, extns) 
