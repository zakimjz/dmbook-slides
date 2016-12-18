#!/usr/bin/env Rscript
library(MASS)
library(e1071)

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
  n = dim(D)[1]
  d = dim(D)[2]
  sampsz = as.integer(n*0.75)
  sampsz = n
  meanL = c(rep(0,3))
  meanQ = c(rep(0,4))
  predsL = c()
  predsQ = c()
  errsL = c(rep(0,folds))
  errsQ = c(rep(0,folds))
  for (f in 1:folds){
	idx = sample(seq(1:n), sampsz, replace=TRUE)
	Z = D[sort(idx),]
    X = Z[,1:(d-1)]
    Y = Z[,d]

	C=10
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
	  AY = c(rep(0,sampsz)) 
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

	  pred <- as.vector(as.integer(predict(model,D[,1:(d-1)])))
	  pred[which(pred==1)]=-1
	  pred[which(pred==2)]=1
	  pp <- pred*D[,d]
	  errors <- which(pp < 0)
	  if (kernel == 'linear'){
		predsL = cbind(predsL,pred)
		errsL[f] = length(errors)/n
	  } else if (kernel == 'quadratic'){
		predsQ = cbind(predsQ,pred)
		errsQ[f] = length(errors)/n
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
  print(errsQ)
  cat(c("mean err L", mean(errsL), " ", sqrt(var(errsL)), "\n"))
  cat(c("mean err Q", mean(errsQ), " ", sqrt(var(errsQ)), "\n"))

  pred = apply(predsL, 1, mode)
  pred[which(pred==2)]=-1
  pp <- pred*D[,d]
  errors <- which(pp < 0)
  cat("error L", length(errors)/n, "\n")

  pred = apply(predsQ, 1, mode)
  pred[which(pred==2)]=-1
  pp <- pred*D[,d]
  errors <- which(pp < 0)
  cat("error Q", length(errors)/n, "\n")

}


############MAIN##########
setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/eval/figs')
set.seed(4100000)
D = read.table('iris-PC.txt')
runs = 30
run_svm(D,runs)

