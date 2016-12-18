#!/usr/bin/env Rscript
library(MASS)
library(mvtnorm)

#D is data, D[,clatr] is labels
run_bootstrap  <- function(D, K){
  n = dim(D)[1]
  clatr = 5
  cls = c("Iris-setosa", "Iris-versicolor", "Iris-virginica")
  err_ary = c(rep(0,K))
  for (f in 1:K){
	idx = sample(seq(1:n), replace=TRUE)
	Tr = D[sort(idx),]
	Te = D
	nTe = dim(Te)[1]
	nTr = dim(Tr)[1]
	D1 = subset(Tr, Tr[,clatr] == "Iris-setosa")
	D2 = subset(Tr, Tr[,clatr] == "Iris-versicolor")
	D3 = subset(Tr, Tr[,clatr] == "Iris-virginica")
	n1 = dim(D1)[1]
	n2 = dim(D2)[1]
	n3 = dim(D3)[1]
  
	cat(n1, n2, n3, "\n")
	#mean and covar matrix
	m1 = colSums(D1[,1:2])/n1
	m2 = colSums(D2[,1:2])/n2
	m3 = colSums(D3[,1:2])/n3
	S1 = (n1-1)/n1*cov(D1[,1:2])
	S2 = (n2-1)/n2*cov(D2[,1:2])
	S3 = (n3-1)/n3*cov(D3[,1:2])

	acc = 0
	err = 0
	Tx = Te
	nx = dim(Tx)[1]
	cl = c(rep(0,nx))
	for (i in 1:nx){
	  x = Tx[i,1:2]
	  p1 = n1/nTr*dmvnorm(x,mean=m1,sigma=S1)
	  p2 = n2/nTr*dmvnorm(x,mean=m2,sigma=S2)
	  p3 = n3/nTr*dmvnorm(x,mean=m3,sigma=S3)
	  idx = which.max(c(p1,p2,p3))
	  cl[i] = cls[idx]
	  if (cl[i] == Tx[i,clatr]){
		  acc = acc+1
	  } else {
		err = err + 1
	  }
	}
	cat(f, acc, acc/nx, "\n")
	cat(f, err, err/nx, "\n")
	err_ary[f] = err/nx
  }
  cat(sort(err_ary), "\n")
  mu = mean(err_ary)
  vv = var(err_ary)
  return (list(mu=mu, var=vv))
}


####MAIN######

setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/eval/figs')
#args <- commandArgs()
#print(args)
#seed = as.integer(args[6])
#set.seed(seed)
#set.seed(4100000)
set.seed(1)
K=50
X = read.table('iris.txt', sep=",")
E = run_bootstrap(X, K)
cat(E$mu, E$var, "\n")
