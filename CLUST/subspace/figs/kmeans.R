#!/usr/bin/env Rscript
setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/subspace/figs')

run_kmeans <- function(){
	eps = 0.001
	means = list()
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
			means[[1]] = FM[[1]]
			means[[2]] = FM[[2]]
			means[[3]] = FM[[3]]
		}
	}
	print(means)

	niters = 0
	notconverged = TRUE
	while(notconverged){
		niters = niters + 1
		oldmu = means
	
		#E_step
		W = matrix(0,n,K)
		M = matrix(unlist(means),d,K)
		dd = matrix(0,n,K)
		for (i in 1:K){
			dd[,i] = colSums((t(D)-M[,i])^2)
		}
		for (j in 1:n){
			ci = which.min(dd[j,])
			W[j,ci] = 1
		}
		
		#M_step
		means = list()
		for (i in 1:K){
			mi = colSums(D*W[,i])/sum(W[,i])
			means[[i]] <- mi
		}

		sumd = 0
		notconverged = TRUE
		for (i in 1:K){
			di = oldmu[[i]]-means[[i]]
			diff = di * di #magnitude of diff
			sumd <- sumd + sum(diff)
		}

		print(c("ITER", niters, sumd))
		if (sumd <= eps){
			notconverged = FALSE
		}
	}
	
	M = matrix(unlist(means),d,K)
	dd = matrix(0,n,K)
	for (i in 1:K){
		dd[,i] = colSums((t(D)-M[,i])^2)
	}
	AssignVec = apply(dd,1,which.min)
	return (AssignVec)
}


#########MAIN###########
#fname = "iris-PC2.txt"
fname = "test1.pca.dat"
fname = "test1.dat"

X = read.table(fname)
if (fname == "test1.dat" || fname == "test1.pca.dat"){
  TL <- c(rep(1,100),rep(2,100),rep(3,100))
  X = cbind(X,TL) 
}
d = ncol(X)-1
K = 3
D = as.matrix(X[,1:d])
n = nrow(D)

#best centers found so far
use_fixed = FALSE
FM = list()
FM[[1]] = c(0.3301256, 1.1180706, 4.6055505)
FM[[2]] = c(7.8478957, 0.7246735, 4.2788653)
FM[[3]] = c(1.567903, 1.440511, 5.204630)


AssignVec = run_kmeans()

#check how many match
TL <- X[,(d+1)]
TP <- cbind(TL, AssignVec)
CT <- table(TP[,1],TP[,2])

print(CT)
purity <- 0
for (i in 1:ncol(CT)){
        purity <- purity + max(CT[,i])
}
print(c(purity,n))
purity <- purity/n
cat('Purity', purity, '\n')
