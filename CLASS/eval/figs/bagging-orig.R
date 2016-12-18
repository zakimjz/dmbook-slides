#!/usr/bin/env Rscript
library(MASS)
library(e1071)
library(class)
library(rpart)
library(RWeka)

IDX = c(2,3,4,5,6,9,12,14,15,16,17,18,19,21,22,23,24,25,26,27,
		29,30,32,34,35,36,37,38,39,40,41,42,43,45,46,47,48,51,
		52,53,54,56,57,58,59,60,61,62,63,64,66,67,69,70,75,76,
		78,85,87,88,89,90,91,92,94,95,96,98,100,101,102,103,
		104,107,109,112,113,114,119,120,121,122,124,125,128,
		130,131,132,133,135,136,137,140,141,142,146,147,150)

mode  <- function(x){
  return (as.numeric(names(table(x)[which.max(table(x))])))
}
#homogeneous quadratic phi mapping
phi2h <- function(x,d){
r = c(rep(0,d*(d+1)/2));
k=1;
for (i in 1:d){
	for (j in i:d){
		if (i == j){
			r[k] = x[i]*x[j];
		}
		else{
			r[k] = sqrt(2)*x[i]*x[j];
		}
		k = k+1;
	}
}
return (r)
}


#D is data, D[,clatr] is labels
run_svm  <- function(D, folds){
  d = dim(D)[2]
  fullidx = seq(1:150)
  trainidx = sample(fullidx)[1:50]
  #trainidx = IDX
  #trainidx = unique(sample(fullidx, replace=TRUE))
  print(trainidx)
  testidx = fullidx[!(fullidx %in% trainidx)]
  print(testidx)
  Tr = D[trainidx,]
  n = dim(Tr)[1]
  Te = D[testidx,]
  nTe = dim(Te)[1]
  meanL = c(rep(0,d))
  meanQ = c(rep(0,(d+1)))
  predsL = c()
  predsQ = c()
  errsL = c(rep(0,folds))
  errsQ = c(rep(0,folds))
  fidx = c()
  for (f in 1:folds){
	idx = sample(seq(1:n), n, replace=TRUE)
	Z = Tr[sort(idx),]
    X = Z[,1:(d-1)]
    Y = Z[,d]

	C=1
	for (kernel in c('linear', 'quadratic')){
	  if (kernel == 'linear'){
		model <- svm(X, Y, type="C-classification",
					 kernel="linear", cost=C, scale=F);
	  } else if (kernel == 'quadratic'){
		model <- svm(X, Y, type="C-classification",
					 kernel="polynomial", degree=2, cost=C, gamma=1, 
					 coef0 = 0, scale=F)
	  }
	  #print(model$index)
	  #print(model$coefs)
	  AY = c(rep(0,n)) 
	  AY[model$index] = model$coefs
	  if (kernel == 'linear'){
			  w <- t(X) %*% AY #compute the weight vector
	  } else if (kernel == 'quadratic'){
		  phiX = apply(X, 1, phi2h, 2)
		  w = phiX %*% AY
	  }

	  #w <- t(model$SV) %*% model$coefs
	  b <- model$rho

	  if (kernel == 'linear'){
		str = sprintf("run %d %s %0.2f*x+%0.2f*y+%0.2f\n", f, kernel, 
					  -1*w[1], -1*w[2], b)
		meanL = meanL + c(w,b)
	  } else if (kernel == 'quadratic'){
		str = sprintf("run %d %s %0.2f*x^2+%0.2f*x*y+%0.2f*y^2+%0.2f\n", f, kernel, 
					  -1*w[1], -1*w[2], -1*w[3], b)
		meanQ = meanQ + c(w,b)
	  }
	  cat(str)

	  ZZ = cbind(X,as.character(Y))
	  names(ZZ)[3] = "V3"
	  #fit = DecisionStump(V3 ~ V1+V2, data=ZZ)
	  #print(ZZ)
	  #fit = J48(V3 ~ V1+V2, data=ZZ, 
	#			control = Weka_control(U=TRUE))
	  #pred <- as.vector(as.integer(predict(fit,Te[,1:(d-1)],
	#							   type='class')))
	  #fit = rpart(V3 ~ V1+V2, method='class', data=Z, 
	  #		  control=rpart.control(minsplit=2))
	  #pred <- as.vector(as.integer(predict(fit,Te[,1:(d-1)],
	#									   type='vector')))
	  pred <- as.vector(as.integer(predict(model,Te[,1:(d-1)])))
	  #pred  <- as.vector(as.integer(knn(X, Te[,1:(d-1)], Y)))
	  #print(pred)
	  pred[which(pred==1)]=-1
	  pred[which(pred==2)]=1
	  pp <- pred*Te[,d]
	  errors <- which(pp < 0)
	  if (kernel == 'linear'){
		predsL = cbind(predsL,pred)
		errsL[f] = length(errors)/nTe
	  } else if (kernel == 'quadratic'){
		predsQ = cbind(predsQ,pred)
		errsQ[f] = length(errors)/nTe
	  }
	}
  }
  meanL = meanL/folds
  str = sprintf("linear mean %0.2f*x+%0.2f*y+%0.2f\n", 
				  -1*meanL[1], -1*meanL[2], meanL[3])
  cat(str)
  meanQ = meanQ/folds
  str = sprintf("quadratic mean %0.2f*x^2+%0.2f*x*y+%0.2f*y^2+%0.2f\n",
				  -1*meanQ[1], -1*meanQ[2], -1*meanQ[3], meanQ[4])
  cat(str)

  print(errsL)
  cat(c("mean err L", mean(errsL), " ", var(errsL), "\n"))
  print(errsQ)
  cat(c("mean err Q", mean(errsQ), " ", var(errsQ), "\n"))

  pred = apply(predsL, 1, mode)
  pred[which(pred==2)]=-1
  pp <- pred*Te[,d]
  errors <- which(pp < 0)
  cat("error L", length(errors)/nTe, "\n")

  pred = apply(predsQ, 1, mode)
  pred[which(pred==2)]=-1
  pp <- pred*Te[,d]
  errors <- which(pp < 0)
  cat("error Q", length(errors)/nTe, "\n")

}


############MAIN##########
setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/eval/figs')
args <- commandArgs(TRUE)
runs = as.integer(args[1])
seed <- as.integer(args[2])
#good seed: 110 for SVM
#set.seed(4100000)
set.seed(seed)
D = read.table('iris-PC.txt')
#runs = 10
run_svm(D,runs)

