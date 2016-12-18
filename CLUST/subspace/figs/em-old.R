#!/usr/bin/Rscript 
#setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/subspace/figs')
library(methods)
library(mvtnorm)
library(foreign)
library(MASS)

apply_dmvnorm <- function(X, mu, S, k){
	EE = eigen(S, symmetric=TRUE, only.values=TRUE)
	detS = prod(EE$values)
	print(c("detS", detS, k))
	if (detS < detThreshold){
		#part1 = (2*pi)**(d/2.0)
		#if (detS == 0){
		#	detS = .Machine$double.xmin
		#} else if (detS == Inf){
		#	detS = .Machine$double.xmax
		#}
		#print(c("detSSSS", detS, k))
	    #part2 = detS**0.5
		#dev = t(X)-mu
		#icov = ginv(S)
		#prodd = apply(dev, 2, function(x){ t(x) %*% icov %*% x})
		#part3 = exp(-0.5*prodd)
		#res =part3/(part1*part2)
		S = S + diag(ridgethresh,d,d)
		print(c("newS", det(S)))
		res = dmvnorm(X,mu,S)
	}
	else{
		res = dmvnorm(X,mu,S)
	}
	return (res)
	#EE = eigen(ss, symmetric=TRUE, only.values=TRUE)
	#vv = prod(EE$values)
	#if (vv < detThreshold){
	#	print(c("vv",i,vv))
	#	ss = ss + diag(ridgethresh,d,d)
	#}
}

###Expectation-Step
E_step <- function(X, Theta, Kmeans, Assign=FALSE){
	W <- matrix(0, nrow=n, ncol=K)
	if (Kmeans){
		M = matrix(unlist(Theta$mu),d,K)
		dd = matrix(0,n,K)
		for (i in 1:K){
			dd[,i] = colSums((t(X)-M[,i])^2)
		}
		for (j in 1:n){
			ci = which.min(dd[j,])
			W[j,ci] = 1
		}
	} else{
		for (i in 1:K){
			m = Theta$mu[[i]]
			ss = Theta$sigma[[i]]
			if (d == 1){
				ss= as.matrix(Theta$sigma[[i]])
			}
			dmv = apply_dmvnorm(X, m, ss, i)
			W[,i] = dmv*Theta$prior[i]
		}
		W = W/rowSums(W)
	}
	if (Assign){
		for (j in 1:n){
			#Assign point to cluster with max prob in final step
			ci = which.max(W[j,])
			Clusters[j] <<- ci
		}
	}
	return (W)
}

###Maximization-Step
M_step <- function(X, W, Kmeans){
	means = list()
	vars = list()
	priors <- colSums(W)/n #new prior
	if (Kmeans){
		for (i in 1:K){
			mi = colSums(X*W[,i])/sum(W[,i])
			means[[i]] <- mi
		}
	}
	else{
		for (i in 1:K){
			mi = colSums(X*W[,i])/sum(W[,i]) #new mean
			cX = t(X)-mi
			S <- cX%*% (t(cX)*W[,i])
			vi = S/sum(W[,i]) #new variance
			#vi = vi + diag(c(rep(0.01,d))) #ridge
			means[[i]] <- mi
			if (length(vi) == 1){
				vars[[i]] = vi
			}
			else{
				if (fullsigma){ vars[[i]] = vi }
				else{ vars[[i]] <- diag(diag(vi)) }
			}
		}
	}
	
	return (list(mu=means, sigma=vars, prior=priors))
}


#Check difference of means for convergence
test_convergence <- function(oldmu, mu){
	sumd = 0
	notconverged = TRUE
	for (i in 1:K){
		di = oldmu[[i]]-mu[[i]]
		diff = di * di #magnitude of diff
		sumd <- sumd + sum(diff)
	}
	return (sumd)
}

#Purity of the final clustering
compute_purity <-function(T, P){
	summ = 0
	for (i in 1:K){
		isect = c(rep(0,K))
		for (j in 1:K){
			isect[j] = length(intersect(P[[i]], T[[j]]))
		}
		summ = summ + max(isect)
	}
	return(summ/nrow(D))
}

#MAIN
args <- commandArgs(TRUE)
#datafile = "iris-PC2.txt"
#datafile = "test1.dat"
#datafile = "clusterdata.txt"
datafile = args[1]
K = as.numeric(args[2])

eps = 0.0001 #stopping threshold
if (length(args) >= 3){
	eps = as.numeric(args[3])
}
ridgethresh = 0.0001
if (length(args) >= 4){
	ridgethresh = as.numeric(args[4])
}
detThreshold = 0.000000001
fullsigma=TRUE
Kmeans=FALSE

print(c("ridge", ridgethresh))
scale = 0.01
if (grep("arff", datafile) > 0){
	X = read.arff(datafile)
} else {
	X = read.table(datafile, as.is=TRUE)
}
n = nrow(X)
print(n)
if (datafile == "test1.dat" || datafile == "test1.pca.dat"){
	TL = c(rep(0,n))
	np = 100
	for (i in 1:K){
		TL[((i-1)*np+1):(i*np)] <- c(rep(i,np))
	}
	X = cbind(X,TL) 
}
d = ncol(X)-1
D = cbind(X[,1:d]*scale, X[,(d+1)])


#best centers found so far
use_fixed = FALSE
FM = list()
FM[[1]] = c(0.379171203,  0.056546122, -0.003323449)
FM[[2]] = c(0.53988846, 0.29295599, 0.03055884)
FM[[3]] = c(0.34294427, 0.27967422, 0.05852964)

means = list()
sigmas = list()
priors = c(rep(1/K,K))
notdone = TRUE
siter = 1
while(notdone){
	print(c("siter", siter))

	for (i in 1:K){
		mu = c()
		for (j in 1:d){
			r  =range(D[,j])
			mu = c(mu, runif(1,r[1], r[2]))
		}
		means[[i]] = mu
	}
	if (use_fixed){
		for (i in 1:K){
			means[[i]] = FM[[i]]
		}
	}

	M = matrix(unlist(means),d,K)
	dd = matrix(0,n,K)
	for (i in 1:K){
		dd[,i] = colSums((t(D[,1:d])-M[,i])^2)
	}
	L  <-  apply(dd,1,which.min)
	CL = c(rep(0,K))
	for (i in 1:K){
		CL[i] = length(which(L==i))
	}
	print(CL)
	minsz = min(CL)
	print(c("minsz", minsz, n/20.0))
	if (minsz < 5){
		notdone = TRUE
		use_fixed = FALSE
	} else {
		notdone = FALSE
	}
	siter = siter+1
}
print(means)


Theta = list(mu=means,sigma=sigmas,prior=priors)
W = E_step(as.matrix(D[,1:d]), Theta, TRUE)
Theta = M_step(as.matrix(D[,1:d]), W, Kmeans)

#start EM iterations
notconverged=TRUE
niters = 0
while(notconverged){
	niters = niters + 1
	oldmu = Theta$mu
	
	W = E_step(as.matrix(D[,1:d]), Theta, Kmeans)
	
	Theta = M_step(as.matrix(D[,1:d]), W, Kmeans)
	
	sumd = test_convergence (oldmu, Theta$mu)
	print(c("ITER", niters, sumd))
	if (sumd <= eps){
		notconverged = FALSE
	}
}

print ('Cluster Parameters')
print(Theta)
print(c('Number of Iterations', niters))

#get final clusters
Clusters <- c(rep(0,n))
W = E_step(D[,1:d], Theta, Kmeans, TRUE)

#check how many match
CT <- table(D[,(d+1)],Clusters)

print(CT)
purity <- 0
for (i in 1:ncol(CT)){
        purity <- purity + max(CT[,i])
}

#diagonal purity, max matching
dPurity = 0.0
nCT = nrow(CT)
y = sort.int(CT, decreasing=TRUE, index.return=TRUE)
rr = as.integer((y$ix-1)%%nCT)+1
cc = as.integer((y$ix-1)/nCT)+1
Tr = c(rep(TRUE, nCT))
Pr = c(rep(TRUE, nCT))
for (i in 1:length(y$x)){
	if (Tr[rr[i]] && Pr[cc[i]]){
		dPurity = dPurity + y$x[i]
		Tr[rr[i]] = FALSE
		Pr[cc[i]] = FALSE
	}
}
print(c(purity, dPurity, n))
purity <- purity/n
dPurity <- dPurity/n
cat('Purity', purity, '\n')
cat ('dPurity', dPurity, '\n')
