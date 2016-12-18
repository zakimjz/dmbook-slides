library(rgl)
setwd("/Users/zaki/research/DataMiningBook/dm08/CLUST/density/figs")
X <- read.table('iris-h0.2.surface.obj', nrows=8364)
D <- as.matrix(X[,2:4])
A = unique(D[,1])
B = unique(D[,2])
C = matrix(D[,3], length(A), length(B))
#surface3d(B,A,C, color='white')
surface3d(B*10,A*10,C*10, color='white')
m = mean(D[,3])
F = C
F[F < m] = m
#surface3d(B,A,F, color='white')
surface3d(B*10,A*10,F*10, color='white')
