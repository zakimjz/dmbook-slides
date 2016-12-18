#!/usr/bin/Rscript

library(fpc)
library(validator)
library(clv)
library(clusterSim)

setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/eval/figs')

getext <- function(datafile){
	D <- read.table(datafile, sep=",")
	X <- D[,1:2]
	C <- D[,4]
	T <- D[,3]

	r = extVal(T, C, index="all")
	print(r)

	cl = cls.attrib(X,C)
	#dummy data to create kmeans object, and then replace info with X,C
	x <- mlbench.2dnormals(150,2)
	km <- kmeans(x$x, 3)
	#now replace
	km$cluster = C
	km$centers = cl$cluster.center
	km$size = cl$cluster.size
	km$totss = km$withinss = km$tot.withinss = km$betweenss = NULL
	clinfo = as.kcca(km, X, simple=TRUE)
	r = intVal(clinfo,X)
	print(r)
	db = index.DB(X,C)
	print(c("db", db$DB))

	dd = dist(X)
	fpc = cluster.stats(dd,C)
	print(c("dunn,ent", fpc$dunn, fpc$entropy))
	#print(c("FPC", fpc))
}

print("bad")
datafile = 'irisPCkmeans-bad.clusters.txt'
getext(datafile)

print("good")
datafile = 'irisPCkmeans-good.clusters.txt'
getext(datafile)

