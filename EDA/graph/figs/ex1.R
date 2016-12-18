setwd('/Users/zaki/research/DataMiningBook/dm08/EDA/graph/figs')
library(igraph)
library(Cairo)

g = graph.empty(n=8, directed=F)

g = add.edges(g, c(0,1,0,2,0,3,0,4,1,4,1,5,2,3,3,4,3,6,4,7,6,7))

diameter(g)

#degree centrality
dc = degree(g)
dc

D = shortest.paths(g)
D

#eccentricity centrality
ec = apply(D, 2, max)
ec
1/ec

#closeness centrality
cc = apply(D,2,sum)
cc
1/cc

#betweenness centrality
bc = betweenness(g)
bc
