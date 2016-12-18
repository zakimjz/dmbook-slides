#!/usr/bin/env Rscript
# gradient descent learning for SVMs
library(methods)
library(MASS)
setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/svm/figs')

normsq <- function(x){
	return (sum(x*x))
}


Jw <- function(w){
	w2 = 0.5*normsq(w)
	wxy <- (w %*% t(XY))
	ind <- sapply(wxy, function(v) return (ifelse(v < 1, 1, 0)))
	sum_yx <- sum((1-wxy)*ind)
	return(w2+C*sum_yx)
}

linesearch  <- function(g){
	rhs = normsq(g)/C
	gxy <- (g %*% t(XY))
	wxy <- (w %*% t(XY))
	indu = which(gxy < 0)
	s=0.001
	l=s
	repeat{
		indw = which(wxy < (1+l*gxy))
		ind = intersect(indu,indw)
		sum_yx = as.vector(c(rep(0,d)))
		nsv = length(ind)
		if (nsv > 1){
			sum_yx <- colSums(XY[ind,])
		} else if (nsv == 1){
			sum_yx = XY[ind,]
		}
		val = t(g)%*%sum_yx
		#print(c("ls", l, nsv, ":", val))
		if (nsv == 0 | val >= rhs){
			l = l-s
			break;
		}
		l = l+s;
	}
	return(l)
}

gradient <- function(w){
	wxy <- (w %*% t(XY))
	#show(wxy)
	#ind <- sapply(wxy, function(v) return (ifelse(v < 1, 1, 0)))
	#print(ind)
	ind = which(wxy<1)
	sum_yx = 0
	if (length(ind) > 1){
		sum_yx <- colSums(XY[ind,])
	} else if (length(ind) == 1){
		sum_yx = XY[ind,]
	}
	g= w - C*sum_yx
	return (g)
}

hessian  <- function(w){
	wxy <- (w %*% t(XY))
	ind <- sapply(wxy, function(v) return (ifelse(v < 1, 1, 0)))
	H <- diag(d) + t(X*ind)%*%(X*ind)
	return (as.matrix(H))
}

gradient2  <- function(w){
	wxy <- (w %*% t(XY))
	ind <- sapply(wxy, function(v) return (ifelse(v < 1, 1, 0)))
	S <- t(X*ind)%*%(X*ind)
	sum_yx <- colSums(XY*ind)
	g= w + (2*C*S)%*%w - 2*C*sum_yx
	show(g)
	return (g)
}

#########MAIN###########
args <- commandArgs(TRUE)
print(args)
trainfile = args[1]
testfile = args[2]
C <- as.real(args[3]) #regularization constant

D <- as.matrix(read.table(trainfile))
n <- nrow(D)
d <- ncol(D)

#t0 <- n*C
t0 <- C

X <- D[,1:(d-1)]

Y <- D[,d]
tau = 1
X <- as.matrix(cbind(X,tau))
XY <- X*Y

w <- c(rep(sqrt(C/d),(d)))
print (w)
print(sqrt(normsq(w)))

#start iterative stochastic gradient descent
iters = 0
eps= 0.001

repeat{
	rate <- 1/(iters+t0)
	w = as.vector(w)
	g <- gradient(w)
	#l <- linesearch(g)
	#if (l > 0){rate=l}
	#H <- hessian(w)
	#wn = as.vector(w - rate*ginv(H)%*%g)
	wn <- w - rate*g
	err=  sqrt(normsq(wn-w))
	show(c("err", err, rate))
	if (err < eps){
        break
	}
	cat('iter', iters, "\n")
    iters <- iters + 1
	w = wn
}

cat(c("\n", "Number of Iterations", iters, "\n"))
cat("\nWeight Vector\n")
print(w)

cat("Norm of w", sqrt(normsq(w)), "\n")

#test the model
Test <- as.matrix(read.table(testfile))
ntest <- nrow(Test)
Xt <- Test[,1:(d-1)]
Yt <- Test[,d]
ones <- c(rep(1,ntest))
Xt <- as.matrix(cbind(Xt,ones))


Yh <- t(w) %*% t(Xt)
preds <- sapply(Yh, function(x) return (ifelse(x>0, 1, -1)))
pp <- preds*Yt
errors <- which(pp < 0)
#print(errors)
cat("\nError Rate\n")
print(length(errors)/ntest)

#margin
wxy <- (w %*% t(XY))
print(which(wxy<1))
print(which(wxy>0.99 & wxy < 1.01))
#print(wxy)
