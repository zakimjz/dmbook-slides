#!/usr/bin/env Rscript

setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/subspace/figs')
X <- read.table("iris.txt", sep=",")
pr <- prcomp(X[,1:4])
n = nrow(X)

D <- pr$x[,1:2]

Xa = which(X[,5] == "Iris-setosa")
Xb = which(X[,5] == "Iris-versicolor")
Xc = which(X[,5] == "Iris-virginica")
#add irrelevant dimension
r  <- range(D)
dda <- runif(length(Xa), r[1], r[2])
ddb <- runif(length(Xb), r[1], r[2])
ddc <- runif(length(Xc), r[1], r[2])
Z = matrix(0,n,3)
Z[Xa,] <- cbind(D[Xa,],dda) #dims 1,2
Z[Xb,] <- cbind(D[Xb,1],ddb, D[Xb,2]) #dims 1,3
Z[Xc,] <- cbind(ddc,D[Xc,]) #dims 2,3
#add labels
l <- c(rep(1,150))
l[Xa] = 0
l[Xb] = 1
l[Xc] = 2
Z = cbind(Z,l)
write.table(Z, "iris-PC2.txt", sep=" ",row.names=F,col.names=F)
