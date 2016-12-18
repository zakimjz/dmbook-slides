#!/usr/bin/env Rscript
library(MASS)
library(mvtnorm)

#D is data, D[,clatr] is labels
run_fold  <- function(D, fold=5){
  n = dim(D)[1]
  clatr = 5
  cls = c("Iris-setosa", "Iris-versicolor", "Iris-virginica")

  permut = sample(seq(1:n))
  sampsz = n/fold
  cat(permut, sampsz, "\n")

  Te = D[permut[1:sampsz],]
  T1 = subset(Te, Te[,clatr] == "Iris-setosa")
  T2 = subset(Te, Te[,clatr] == "Iris-versicolor")
  T3 = subset(Te, Te[,clatr] == "Iris-virginica")
  write.table(T1[,1:2], "/tmp/ff.txt",  row.names=FALSE, col.names=FALSE)
  write.table(T2[,1:2], "/tmp/ff.txt",  row.names=FALSE, col.names=FALSE,
                          append=TRUE)
  write.table(T3[,1:2], "/tmp/ff.txt",  row.names=FALSE, col.names=FALSE,
                          append=TRUE)
  cat(dim(T1), dim(T2), dim(T3), "\n")



  Tr = D[permut[(sampsz+1): n],]
  nTe = dim(Te)[1]
  nTr = dim(Tr)[1]
  cat("TrTe", nTe, nTr, "\n")
  flush.console()
  #print(Te[,clatr])
  D1 = subset(Tr, Tr[,clatr] == "Iris-setosa")
  D2 = subset(Tr, Tr[,clatr] == "Iris-versicolor")
  D3 = subset(Tr, Tr[,clatr] == "Iris-virginica")
  n1 = dim(D1)[1]
  n2 = dim(D2)[1]
  n3 = dim(D3)[1]
 
  cat("size", n1, n2, n3, "\n")
  #mean and covar matrix
  m1 = colSums(D1[,1:2])/n1
  m2 = colSums(D2[,1:2])/n2
  m3 = colSums(D3[,1:2])/n3
  S1 = (n1-1)/n1*cov(D1[,1:2])
  S2 = (n2-1)/n2*cov(D2[,1:2])
  S3 = (n3-1)/n3*cov(D3[,1:2])
  #S1 = diag(diag(S1))
  #S2 = diag(diag(S2))
  #S3 = diag(diag(S3))
  print(m1)
  print(m2)
  print(m3)
  print(S1)
  print(S2)
  print(S3)
  print("inverses")
  print(ginv(S1))
  print(ginv(S2))
  print(ginv(S3))

  acc = 0
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
	}
  }
  cat(acc, acc/nx, "\n")
  CT = table(cl, Tx[,clatr])
  print("contingency table")
  print(CT)
}


####MAIN######

setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/eval/figs')
#args <- commandArgs()
#print(args)
#seed = as.integer(args[6])
#set.seed(seed)
#set.seed(4100000)
set.seed(1)

X = read.table('iris.txt', sep=",")
run_fold(X)
Z = cbind(X[,1:2], as.integer(X[,5]))
D1 = subset(Z, Z[,3] == 1)
D2 = subset(Z, Z[,3] == 2)
D3 = subset(Z, Z[,3] == 3)

write.table(D1, "/tmp/ffd.txt",  row.names=FALSE, col.names=FALSE)
write.table(D2, "/tmp/ffd.txt",  row.names=FALSE, col.names=FALSE,
			append=TRUE)
write.table(D3, "/tmp/ffd.txt",  row.names=FALSE, col.names=FALSE,
			append=TRUE)
