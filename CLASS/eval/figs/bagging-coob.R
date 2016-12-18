#!/usr/bin/env Rscript
library(MASS)
library(e1071)
library(class)
library(rpart)
library(RWeka)
library(ipred)
#library(adabag)

stump=rpart.control(cp=-1,maxdepth=1,minsplit=0)
default=rpart.control()
pow2=rpart.control(cp=-1,maxdepth=30,minsplit=0)

mymode  <- function(x){
  #return (as.numeric(names(table(x)[which.max(table(x))])))
  return(sign(sum(x)))
}

numzeros <- function(x){
  return(length(which(x == 0)))
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
  n = dim(D)[1]
  d = dim(D)[2]
  print(dim(D))
  meanL = c(rep(0,d))
  meanQ = c(rep(0,(d+1)))
  predsL = c()
  predsQ = c()
  errsL = c(rep(0,folds))
  errsQ = c(rep(0,folds))
  fullidx = seq(1:n)
  for (f in 1:folds){
	idx = sample(seq(1:n), n, replace=TRUE)
    testidx = fullidx[!(fullidx %in% idx)]
	Te = D[testidx,]
	nTe = dim(Te)[1]

	Z = D[sort(idx),]
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
	  fit = rpart(V3 ~ V1+V2, method='class', data=ZZ, 
	  		  control=pow2)
	  pred <- as.vector(as.integer(predict(fit,Te[,1:(d-1)],
	   									   type='vector')))
	  #pred <- as.vector(as.integer(predict(model,Te[,1:(d-1)])))
	  #pred  <- as.vector(as.integer(knn(X, Te[,1:(d-1)], Y)))
	  #print(pred)
	  pred[which(pred==1)]=-1
	  pred[which(pred==2)]=1
	  pp <- pred*Te[,d]
	  errors <- which(pp < 0)

	  predX = c(rep(0,n))
	  for (i in 1:n){
		wi = which(testidx == i)
		if (length(wi) > 0){
		  predX[i] = pred[wi]
		}
	  }
	  
	  if (kernel == 'linear'){
		predsL = cbind(predsL,predX)
		errsL[f] = length(errors)/nTe
	  } else if (kernel == 'quadratic'){
		predsQ = cbind(predsQ,predX)
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

  pred = apply(predsL, 1, mymode)
  ZZZ = apply(predsL,1,numzeros)
  omitidx = which(ZZZ == runs)
  valididx = fullidx[!(fullidx %in% omitidx)]
  pred[which(pred==0)]=1
  pp <- pred[valididx]*D[valididx,d]
  errors <- which(pp < 0)
  cat("error L", length(errors)/length(valididx), "\n")

  pred = apply(predsQ, 1, mymode)
  ZZZ = apply(predsQ,1,numzeros)
  omitidx = which(ZZZ == runs)
  valididx = fullidx[!(fullidx %in% omitidx)]
  pred[which(pred==0)]=1
  pp <- pred[valididx]*D[valididx,d]
  errors <- which(pp < 0)
  cat("error Q", length(errors)/length(valididx), "\n")

  ZZ = cbind(D[,1:(d-1)],as.character(D[,d]))
  names(ZZ)[3] = "V3"

  res = bagging(V3~V1+V2, data=ZZ, nbagg=runs, coob=TRUE,
				control=pow2)
  print(res$err)
  #pred = as.vector(as.integer(res$y))
  #print(pred)
  #pred[which(pred==1)]=-1
  #pred[which(pred==2)]=1
  #pp <- pred*D[,d]
  #print(pp)
  #errors <- which(pp < 0)
  #cat("error bag", length(errors)/n, "\n")

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

