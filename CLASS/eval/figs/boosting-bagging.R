#!/usr/bin/env Rscript
library(MASS)
library(e1071)
library(kernlab)
library(class)
library(rpart)
library(RWeka)
library(ipred)
#library(adabag)


normsq <- function(x){
        return (sum(x*x))
}

gradient <- function(w, XY, Cvar){
        wxy <- (w %*% t(XY))
        ind = which(wxy<1)
        sum_yx = 0
        if (length(ind) > 1){
                sum_yx <- colSums(XY[ind,])
        } else if (length(ind) == 1){
                sum_yx = XY[ind,]
        }
        g= w - Cvar*sum_yx
        return (g)
}

stump=rpart.control(cp=-1,maxdepth=1,minsplit=0)
default=rpart.control()
pow2=rpart.control(cp=-1,maxdepth=30,minsplit=0)

mymode  <- function(x){
  #print(x)
  #tt = table(x)
  #print(tt)
  #print(names(tt))
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
run_boost  <- function(Tridx, Teidx, K, alg, C=1){
  n = length(Tridx)
  #print(Tridx)
  str = ""
  ww = c()
  predsTe = matrix(0,length(Teidx),K)
  predsTr = matrix(0,length(Tridx),K) 
  MxxTr = matrix(0,length(Tridx),K)
  MxxTe = matrix(0,length(Teidx),K)
  errs = c(rep(0,K))
  delTr = c(rep(0,K))
  delTe = c(rep(0,K))
  meanH = c(rep(0,(d+1)))
  wt = c(rep(1/n, n))
  t = 1
  alphat = c(rep(0,K))
  numrestarts = 0
  numbadSVMs = 0
  while(t <= K){
	#weighted sample
	#print(c(t, " --", wt, "\n"))
	idx = sample(Tridx, n, replace=TRUE, prob=wt)
	#print(sort(idx))
    X = Z[idx,1:(d-1)]
    Y = Z[idx,d]
	ZZ = Z[idx,]
	if (alg == 'linearSVM'){
	  model <- svm(X, Y, type="C-classification",
				   kernel="linear", cost=C, scale=FALSE);
	  AY = c(rep(0,n)) 
	  AY[model$index] = model$coefs
	  w <- t(X) %*% AY #compute the weight vector
	  b <- model$rho
	  ww = c(-w[1], -w[2], b, 0)
	  #print(ww)
	  str = sprintf("%0.2f*x+%0.2f*y+%0.2f", 
					-1*w[1], -1*w[2], b)
	  #cat(str)
	  pred <- as.vector(as.integer(predict(model,D[Tridx,1:(d-1)])))
	  predTe <- as.vector(as.integer(predict(model,D[Teidx,1:(d-1)])))
	} else if (alg == 'primalSVM'){
	  X <- D[idx,1:(d-1)]
	  Y <- D[idx,d]
	  tau = 1
	  X <- as.matrix(cbind(X,tau))
	  XX <- as.matrix(cbind(D[,1:(d-1)],tau))
	  XY <- X*Y
	  t0 <- C
	  w <- c(rep(sqrt(C/d),(d)))
	  #start iterative stochastic gradient descent
	  iters = 0
	  eps= 0.001
	  repeat{
		rate <- 1/(iters+t0)
		w = as.vector(w)
		g <- gradient(w, XY, C)
		wn <- w - rate*g
		err=  sqrt(normsq(wn-w))
		#show(c("err", err, rate))
		if (err < eps){
		break
		}
		#cat('iter', iters, "\n")
		iters <- iters + 1
		w = wn
	  }
	  #print(c(t, "--", w))
	  #print(c("www", w))
	  ww = c(w[1], w[2], w[3], 0)
	  b = w[3]
	  #print(ww)
	  str = sprintf("%0.2f*x+%0.2f*y+%0.2f", 
					w[1], w[2], b)
	  #cat(str)
	  Yh <- t(w) %*% t(XX[Tridx,])
	  pred <- sapply(Yh, function(x) return (ifelse(x>0, 1, -1)))
	  YhTe <- t(w) %*% t(XX[Teidx,])
	  predTe <- sapply(YhTe, function(x) return (ifelse(x>0, 1, -1)))
	} else if (alg == 'linearSVMkernlab'){
	  model <- ksvm(V3~V1+V2, data=ZZ, type="C-svc",
				   kernel="vanilladot", C=C, scaled=FALSE)
	  #print(model)
	  AY = c(rep(0,n)) 
	  aidx = alphaindex(model)[[1]]
	  AY[aidx] = coef(model)[[1]]
	  #print(AY)
	  w <- t(X) %*% AY #compute the weight vector
	  b <- b(model)
	  ww = c(-w[1], -w[2], b, 0)
	  #print(ww)
	  str = sprintf("%0.2f*x+%0.2f*y+%0.2f", 
					-1*w[1], -1*w[2], b)
	  #cat(str)
	  pred <- as.vector(as.integer(predict(model,D[Tridx,1:(d-1)])))
	  predTe <- as.vector(as.integer(predict(model,D[Teidx,1:(d-1)])))
	} else if (alg == 'quadraticSVM'){
	  model <- svm(X, Y, type="C-classification",
				   kernel="polynomial", degree=2, cost=C, gamma=1, 
				   coef0 = 0, scale=FALSE)
	  AY = c(rep(0,n)) 
	  AY[model$index] = model$coefs
	  phiX = apply(X, 1, phi2h, 2)
	  w = phiX %*% AY
	  b <- model$rho
	  ww = c(-w[1], -w[2], -w[3], b)
	  str = sprintf("%0.2f*x^2+%0.2f*x*y+%0.2f*y^2+%0.2f", 
					-1*w[1], -1*w[2], -1*w[3], b)
	  #cat(str)
	  pred <- as.vector(as.integer(predict(model,D[Tridx,1:(d-1)])))
	  predTe <- as.vector(as.integer(predict(model,D[Teidx,1:(d-1)])))
	  phiTr <- apply(D[Tridx, 1:(d-1)], 1, phi2h, 2)
	  phiTe <- apply(D[Teidx, 1:(d-1)], 1, phi2h, 2)
	} else if (alg == "DecisionStump"){
	  fit = DecisionStump(V3 ~ V1+V2, data=ZZ)
	  pred <- as.vector(as.integer(predict(fit,D[Tridx,1:(d-1)],
								   type='class')))
	  predTe <- as.vector(as.integer(predict(fit,D[Teidx,1:(d-1)],
											 type="class")))
	} else if (alg == "J48"){
	  fit = J48(V3 ~ V1+V2, data=ZZ, 
	  			control = Weka_control(U=TRUE))
	  pred <- as.vector(as.integer(predict(fit,D[Tridx,1:(d-1)],
								   type='class')))
	  predTe <- as.vector(as.integer(predict(fit,D[Teidx,1:(d-1)],
											 type="class")))
	} else if (alg == "rpartFull"){
	  fit = rpart(V3 ~ V1+V2, method='class', data=ZZ, 
	  		control=pow2)
	  pred <- as.vector(as.integer(predict(fit,D[Tridx,1:(d-1)],
	   									 type='vector')))
	  predTe <- as.vector(as.integer(predict(fit,D[Teidx,1:(d-1)],
											 type="vector")))
	} else if (alg == "rpartStump"){
	  fit = rpart(V3 ~ V1+V2, method='class', data=ZZ, 
	  		control=stump)
	  pred <- as.vector(as.integer(predict(fit,D[Tridx,1:(d-1)],
	   									 type='vector')))
	  predTe <- as.vector(as.integer(predict(fit,D[Teidx,1:(d-1)],
											 type="vector")))
	} else if (alg == "knn"){
	  pred  <- as.vector(as.integer(knn(X, D[Tridx,1:(d-1)], Y)))
	  predTe <- as.vector(as.integer(knn(X,D[Teidx,1:(d-1)],Y)))
	}
	#print(c("pred",pred))
	pred[which(pred==1)]=-1
	pred[which(pred==2)]=1
	predTe[which(predTe==1)]=-1
	predTe[which(predTe==2)]=1
	#print(pred)

	pp <- pred*D[Tridx,d]
	#print(pp)
	Ivec = c(rep(0,n))
	Ivec[which(pp < 0)] = 1 #set misclassified to 1
	#print(c("Ivec", Ivec))
	et  <- wt %*% Ivec
	errs[t] = et
	if (et > 0.5){
	  pred = -1*pred
	  predTe = -1*predTe
	}
	predsTr[,t] = pred
	predsTe[,t] = predTe
	#print(c("error", t, et))
	badSVM = FALSE
	if (length(grep("SVM",alg)) > 0){
	  if (abs(ww[1]) < 0.01 && abs(ww[2]) < 0.01 
		&& abs(1-abs(ww[3])) < 0.01){
		numbadSVMs = numbadSVMs + 1
		badSVM = TRUE
	  } 
	  else{
		meanH = meanH + ww
	  }
	}
	if (bagging){
	  wt = c(rep(1/n, n))
	  alphat[t] = 1
	  t = t+1
	  if (length(grep("SVM",alg)) > 0){
		if (badSVM){
			t = t-1
		} else {
		  strx = sprintf('\\psplotImp[algebraic](-4.1,-3.1)(4.1,3.1){%s}',
					   str)
		  print(strx)
		} 
	  }
	} else {
	  if (et == 0){ 
		break 
	  } else {
		alphat[t] = log((1-et)/et)
		if (et < 0.5){
		  if (length(grep("SVM",alg)) > 0){
			strx = sprintf(
			  '\\psplotImp[algebraic](-4.1,-3.1)(4.1,3.1){%s}',
			  str)
			print(strx)
		  }
		  Uvec = Ivec*(1-et)/et
		  Uvec[which(Uvec == 0)] = 1
		  #print(c("Uvec", Uvec))
		  wt = wt*Uvec
		  wt = wt / sum(wt)
		  t = t+1
		} else {
		  numrestarts = numrestarts + 1
		}
	  }
	}
  }
  nK = t-1 #new K value
  #normalize alphat
  alphat = alphat/sum(alphat)
  print(c("alpha", alphat))
  print(c("errs", errs))
  #print(dim(predsTr))
  #print(predsTr)
  WpredsTr = t(predsTr)*alphat
  wpredTr = apply(WpredsTr, 2, mymode)
  #print(wpred)
  pp <- wpredTr*D[Tridx,d]
  errary = which(pp < 0) #set misclassified to 1
  errTr = length(errary)/length(Tridx)
  sqerrTr = mean((wpredTr - D[Tridx,d])**2)
  print(c("final Tr error", errTr, sqerrTr))
  meanH = meanH/nK
  meanstr = sprintf("%0.3f %0.3f %0.3f %0.3f", meanH[1], meanH[2],
					meanH[3], meanH[4])
  print(c("meanH", nK, meanstr) )
  WpredsTe = t(predsTe)*alphat
  wpredTe = apply(WpredsTe, 2, mymode)
  pp <- wpredTe*D[Teidx,d]
  errary = which(pp < 0) #set misclassified to 1
  errTe = length(errary)/length(Teidx)
  sqerrTe = mean((wpredTe - D[Teidx,d])**2)
  print(c("final Te error", errTe, sqerrTe))
  print(numrestarts)
  print(numbadSVMs)

  #compute squared-loss error, avg bias and avg var
  lTr = (predsTr - D[Tridx,d])**2
  lossTr = mean(lTr)
  mTr = apply(predsTr,1,mean) #mean pred value for each point
  AVTr = apply((predsTr - mTr)**2, 1, mean) #var for each point
  ABTr = (mTr - D[Tridx,d])**2 #bias for each point
  muAVTr = mean(AVTr)
  muABTr = mean(ABTr)
  
  lTe = (predsTe - D[Teidx,d])**2
  lossTe = mean(lTe)
  mTe = apply(predsTe, 1, mean)
  AVTe = apply((predsTe - mTe)**2, 1, mean)
  ABTe = (mTe - D[Teidx,d])**2	
  muAVTe = mean(AVTe)
  muABTe = mean(ABTe)
  return(c(errTr, errTe, sqerrTr, sqerrTe, muAVTr, muABTr, muAVTe, muABTe,
		   lossTr, lossTe))
}

run_CVboost <- function(K,alg,fold=5, C=1){
  fullidx = seq(1:N)
  permut = sample(seq(1:N))
  sampsz = as.integer(N/fold)
  ERRS = matrix(0,(fold+1),10)
  if (usefulldata){
	ff = c(fold+1)
  } else {
	ff = seq((fold+1),1)
  }
  for (f in ff){
	if (f == (fold+1)){ #train and test on full data
	  Teidx = permut
	  Tridx = permut
	} else {
	  lb = (f-1)*sampsz+1
	  ub = f*sampsz
	  if (f==fold){ub = N}
	  Teidx = permut[lb:ub]
	  Tridx = fullidx[!(fullidx %in% Teidx)]
	}
	#print(Teidx)
	#print(Tridx)
	eb = run_boost(Tridx, Teidx, K, alg, C) 
	print(c(f, eb))
	ERRS[f,] = eb
  }
  print(ERRS)
  means = colSums(ERRS[1:fold,])/fold
  return(c(K, means, ERRS[(fold+1),]))
}


run_diff_K  <- function(runs,alg,skip=1){
  RES = c()
  for (K in seq(0,runs,skip)){
	if (K == 0){K = K+1}
	set.seed(seed)
	res = run_CVboost(K,alg)
	RES = rbind(RES,res)
  }
  print(RES)
}

run_diff_C  <- function(K, alg){
  RES = c()
  cvals = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)
  for (Cvar in cvals){
	set.seed(seed)
	res = run_CVboost(K,alg,C=Cvar)
	RES = rbind(RES,res)
  }
  RES = cbind(cvals,RES)
  print(RES)
}

############MAIN##########
setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/eval/figs')
args <- commandArgs(TRUE)
K = as.integer(args[1])
alg <- args[2]
seed <- as.integer(args[3])
bagging <- FALSE
usefulldata=FALSE
use_diff_K= FALSE
if (length(args) > 3){
  bagging <- as.integer(args[4])
}
if (length(args) > 4){
  usefulldata <- as.integer(args[5])
}
if (length(args) > 5){
  use_diff_K = TRUE
  skip <- as.integer(args[6])
}

print(c(K,alg,seed,bagging,usefulldata))
#good seed: 110 for SVM
#set.seed(4100000)
set.seed(seed)
D = read.table('iris-PC.txt')
N = dim(D)[1]
d = dim(D)[2]
Z = cbind(D[,1:(d-1)],as.character(D[,d]))
names(Z)[3] = "V3"
#if (use_diff_K){
#  run_diff_K(K,alg,skip) #comment this out to vary K
#} else {
#  run_CVboost(K,alg)
#}
run_diff_C(K,alg)
