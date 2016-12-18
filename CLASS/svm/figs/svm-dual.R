#!/usr/bin/env Rscript
# gradient descent learning for SVMs
library(methods)
compute_rate <- function(maxNorm){
    R <- c(rep(0,n))
	dot <- diag(K)
	#print(dot)
	R <- 1/dot
	#print(R)
	maxNorm <<- max(dot)
	maxNorm <<- sqrt(maxNorm)

    return(R)
}

normsq <- function(x){
        return (sum(x*x))
}

get_dot_sum <- function(idx){
	AY <- A*Y
	ssum <- sum(AY*K[,idx])
	return (ssum)
}

projectval <- function(a){
	if (a < 0) return (0)
	else if (a > C) return (C)
	else return (a)
}

#homogeneous quadratic phi mapping
phi2h <- function(x){
	#since points are in 1 higher dimension, the individual
	#dimenson linear components are obtained by dot-prod with the
	#last dimension, which is 1
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

kernel = 'linear'
if (length(args) > 3){
	kernel = args[4]
}

D <- as.matrix(read.table(trainfile))
n <- nrow(D)
d <- ncol(D)

X <- D[,1:(d-1)]
Y <- D[,d]
ones <- c(rep(1,nrow(D)))
X <- as.matrix(cbind(X,ones))

if (kernel == 'linear'){
	K <- X %*% t(X);
} else if (kernel == 'quadratic'){
	K <- X %*% t(X);
	#K = 1+K;
	K <- K*K;
}

A <- c(rep(0,n))
maxNorm <- 0
R <- compute_rate(maxNorm)

#start iterative stochastic gradient ascent
iters = 0
eps = 0.00001
repeat{
	#idxs = sample(n)
	Aold = A
	for (k in 1:n){
		#i = idxs[k]
		i = k
		g <- 1-Y[i]*get_dot_sum(i)
		A[i] <- A[i] + R[i]*g
		A[i] <- projectval(A[i])
    }
    iters <- iters + 1
	err = sqrt(normsq(A-Aold))
	show(c("err", err))
    if (err < eps ){
        break
	}
}

cat(c("\n", "Number of Iterations", iters, "\n"))
sv <- sapply(A, function(x) return (ifelse (x>0, 1, 0)))
nsv <- sum(sv)
cat("\nNumber of support vectors", nsv, "\n")
cat("\nSupport Vectors\n")
print(which(A>0))
cat("\nAlpha Values\n")
print(A[which(A>0)])

#test the model
Test <- as.matrix(read.table(testfile))
ntest <- nrow(Test)
Xt <- Test[,1:(d-1)]
Yt <- Test[,d]
ones <- c(rep(1,ntest))
Xt <- as.matrix(cbind(Xt,ones))
Kt <- X %*% t(Xt)

if (kernel == 'linear'){
	Kt <- X %*% t(Xt);
} else if (kernel == 'quadratic'){
	Kt <- X %*% t(Xt);
	#Kt = 1+Kt;
	Kt <- Kt*Kt;
}

Yh <- t(Kt)%*%(A*Y)

preds <- sapply(Yh, function(x) return (ifelse(x>0, 1, -1)))
pp <- preds*Yt
errors <- which(pp < 0)
#print(errors)
cat("\nError Rate\n")
print(length(errors)/ntest)

if (kernel == 'linear'){
	w <- t(X) %*% (A*Y) #compute the weight vector
} else if (kernel == 'quadratic'){
    phiX = apply(X, 1, phi2h)
    w = phiX %*% (A*Y)
}
cat("\nWeight Vector\n")
print(c(w))

#if (kernel == 'linear'){
	#ind = A[which(A>0 & A < C)]
	#wxy <- (X[ind,]*Y[ind]) %*% w
	#print(wxy)
	#margin = min(wxy)
	#print(margin)
	#w = w*(1/margin)
	#print(w)
#}

