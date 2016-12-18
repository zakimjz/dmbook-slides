#!/usr/bin/env Rscript
# gradient descent learning for SVMs
library(methods)
library(MASS)
setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/svm/figs')

normsq <- function(x){
	return (sum(x*x))
}


gradient <- function(b){
	yKb <- b %*% t(K*Y)
	#print(yKb)
	#ind <- sapply(yKb, function(v) return (ifelse(v < 1, 1, 0)))
	#print(ind)
	#sum_yx <- as.vector(colSums(K*Y*ind))
	#show(sum_yx)	
	ind = which(yKb<1)
	sum_yx = c(rep(0,n))
	if (length(ind) > 1){
		sum_yx <- colSums((K*Y)[ind,])
	} else if (length(ind) == 1){
		sum_yx = (K*Y)[ind,]
	}
	#show(sum_yx)	
	#print(Y)
	#print(K)
	#print(sum_yx)
	g= K%*%b - C*sum_yx
	g = as.vector(g)
	#show(g)	
	return (g)
}

projectval <- function(a){
        if (a < -C) return (-C)
        else if (a > C) return (C)
        else return (a)
}



#homogeneous quadratic phi mapping
phi2h <- function(x){
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


#########MAIN###########
args <- commandArgs(TRUE)
print(args)
trainfile = args[1]
testfile = args[2]
C <- as.real(args[3]) #regularization constant

options(digits=3)
kernel = 'linear'
if (length(args) > 3){
    kernel = args[4] #only quadratic or linear
}


D <- as.matrix(read.table(trainfile))
n <- nrow(D)
d <- ncol(D)


X <- D[,1:(d-1)]

Y <- D[,d]
tau = 1
X <- as.matrix(cbind(X,tau))

if (kernel == 'linear'){
        K <- X %*% t(X);
}

if (kernel == 'quadratic'){
        K <- X %*% t(X);
        K <- K*K;
}


KI = ginv(K)

b <- c(rep(sqrt(C/d),n))
print (b)
print(sqrt(normsq(b)))

#start iterative stochastic gradient descent
iters = 0
eps= 0.00001
maxiters = 100000

KI = ginv(K)
mk = max(K)
#t0 <- C
t0 <- C

repeat{
	rate <- 1/(iters+t0)
	b = as.vector(b)
	g <- gradient(b)
	#H <- hessian(w)
	#rate = normsq(g)/(t(g)%*%K%*%g)
	#bn <- b - rate/mk*g
	
	bn = as.vector(b - rate*(KI%*%g))

	err=  sqrt(normsq(bn-b))
	show(c("err",iters, err, rate))
	if (err < eps | iters > maxiters){
        break
	}
    iters <- iters + 1
	b = bn
}

cat(c("\n", "Number of Iterations", iters, "\n"))
cat("\nCoeff Vector\n")
print(b)

cat("Norm of b", sqrt(normsq(b)), "\n")

if (kernel == 'linear'){
	w = t(X) %*% b
}
if (kernel == 'quadratic'){
	phiX = apply(X, 1,phi2h)
	w = phiX %*% b
}

show(c(w))
show(sqrt(norm(w)))

#test the model
Test <- as.matrix(read.table(testfile))
ntest <- nrow(Test)
Xt <- Test[,1:(d-1)]
Yt <- Test[,d]
ones <- c(rep(1,ntest))
Xt <- as.matrix(cbind(Xt,ones))

if (kernel == 'quadratic'){
	Xt = t(apply(Xt,1,phi2h))
}

Yh <- t(w) %*% t(Xt)
preds <- sapply(Yh, function(x) return (ifelse(x>0, 1, -1)))
pp <- preds*Yt
errors <- which(pp < 0)
#print(errors)
cat("\nError Rate\n")
print(length(errors)/ntest)

