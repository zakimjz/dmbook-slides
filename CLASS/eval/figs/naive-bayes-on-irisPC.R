#!/usr/bin/env Rscript
library(MASS)
library(mvtnorm)
library(ROCR)

print_tpr_fpr  <- function(TP, FP, nP, nN){
  tpr = TP/nP
  fpr = FP/nN
  cat(fpr, tpr, "\n")
}

AUC <- function(oldFP, oldTP, FP, TP, nP, nN){
  ofpr = oldFP/nN
  otpr = oldTP/nP
  fpr = FP/nN
  tpr = TP/nP
  return (abs(fpr-ofpr)*(tpr+otpr)/2)
}

roc_plot  <-  function(P, Te, clatr){
  Sidx = sort.int(P, decreasing=TRUE, index.return=TRUE)
  P = P[Sidx$ix]
  Te = Te[Sidx$ix,]
  #cat(format.default(P, digits=3, sci=FALSE), "\n")
  cat(P,"\n")
  print(Te[,clatr])
  TP = FP = 0
  T1 = subset(Te, Te[,clatr] == -1)
  T2 = subset(Te, Te[,clatr] == 1)
  nP = dim(T1)[1]
  nN = dim(T2)[1]
  nx = dim(Te)[1]
  theta = 100
  A = 0
  oldTP = oldFP = 0
  for (i in 1:nx){
	if (theta > P[i]){
	  print_tpr_fpr(TP, FP, nP, nN)
	  A = A+ AUC(oldFP, oldTP, FP, TP, nP, nN)
	  cat("FPTP", i, FP, oldFP, TP, oldTP, "\n")
	  theta = P[i]
	  oldFP = FP
	  oldTP = TP
	}
	if (Te[i,clatr] == -1){ #positive class
		TP = TP + 1
	} else {
		FP = FP + 1
	}
  }
  print_tpr_fpr(TP, FP, nP, nN)
  A = A+ AUC(oldFP, oldTP, FP, TP, nP, nN)
  #A = A/(nP * nN)
  cat("AUC", A, "\n")
}


#D is data, D[,clatr] is labels
run_fold  <- function(D, fold=5){
  n = dim(D)[1]
  clatr = 3

  permut = sample(seq(1:n))
  sampsz = n/fold
  cat(permut, sampsz, "\n")

  Te = D[permut[1:sampsz],]
  T1 = subset(Te, Te[,clatr] == -1)
  T2 = subset(Te, Te[,clatr] == 1)
  write.table(T1, "/tmp/ff.txt",  row.names=FALSE, col.names=FALSE)
  write.table(T2, "/tmp/ff.txt",  row.names=FALSE, col.names=FALSE,
			  append=TRUE)
  cat(dim(T1), dim(T2), "\n")

  Tr = D[permut[(sampsz+1): n],]
  nTe = dim(Te)[1]
  nTr = dim(Tr)[1]
  cat("TrTe", nTe, nTr, "\n")
  #print(Te[,clatr])
  D1 = subset(Tr, Tr[,clatr] == -1)
  D2 = subset(Tr, Tr[,clatr] == 1)
  n1 = dim(D1)[1]
  n2 = dim(D2)[1]
 
  cat("size", n1, n2, "\n")
  #mean and covar matrix
  m1 = colSums(D1[,1:2])/n1
  m2 = colSums(D2[,1:2])/n2
  S1 = (n1-1)/n1*cov(D1[,1:2])
  S2 = (n2-1)/n2*cov(D2[,1:2])
  S1 = diag(diag(S1))
  S2 = diag(diag(S2))
  print(m1)
  print(m2)
  print(S1)
  print(S2)
  print("inverses")
  print(ginv(S1))
  print(ginv(S2))

  acc = 0
  Tx = Te
  nx = dim(Tx)[1]
  cl  = c(rep(0,nx))
  pr  = c(rep(0,nx))
  for (i in 1:nx){
	x = Tx[i,1:2]
	p1 = n1/nTr*dmvnorm(x,mean=m1,sigma=S1)
	p2 = n2/nTr*dmvnorm(x,mean=m2,sigma=S2)
	pr[i] = p1/(p1+p2)
	if (p1 > p2) {
	  cl[i] = -1
	} else {
	  cl[i] = 1
	}
	if (cl[i] == Tx[i,clatr]){
		acc = acc+1
	}
  }
  cat(acc, acc/nx, "\n")
  CT = table(cl, Tx[,clatr])
  print("contingency table")
  print(CT)

  roc_plot(pr, Tx, clatr)

  Sidx = sort.int(pr, decreasing=TRUE, index.return=TRUE)
  #print(cbind(Tx[Sidx$ix, clatr],pr[Sidx$ix]))
  pred <- prediction(pr[Sidx$ix], Tx[Sidx$ix,clatr], label.ordering = c(1,-1))
  #print(pred)
  perf <- performance(pred, measure='tpr', x.measure='fpr')
  perf2 <- performance(pred, measure='auc')

  print(perf2)
  #plot(perf)

  #SVM
  svmf = ksvm(V[,clatr]~.,data=Tx)

}


####MAIN######

setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/eval/figs')
#args <- commandArgs()
#print(args)
#seed = as.integer(args[6])
#set.seed(seed)
set.seed(4100000)
X = read.table('iris-PC.txt')
run_fold(X)
