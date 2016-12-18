#!/usr/bin/Rscript 
setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/representative/figs')
library(methods)
library(mvtnorm)

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
			W[,i] = dmvnorm(X, m, ss)*Theta$prior[i]
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
#datafile = 'data3c.txt'
#datafile = 'iris.txt'
datafile = 'iris-PC.txt'
datafile = 'em-1d-data.txt'
#datafile = 'kerneldata.txt'
#datafile = 'foo.txt'
K =2
fullsigma=TRUE
Kmeans=FALSE
eps = 0.0001 #stopping threshold

D = read.table(datafile, sep=",")
d = ncol(D)-1
n = nrow(D)


means = list()
sigmas = list()
priors = c(rep(1/K,K))
for (i in 1:K){
	mu = c()
	sigma = c()
	for (j in 1:d){
		r  =range(D[,j])
		mu = c(mu, runif(1,r[1], r[2]))
		#sigma = c(sigma, runif(1,r[1], r[2]))
		if (d == 1){
			sigma = 1
		} else {
			sigma = diag(d)
		}
	}
	means[[i]] = mu
	sigmas[[i]] = sigma
}
Theta = list(mu=means,sigma=sigmas,prior=priors)
print(Theta)

#start EM iterations
notconverged=TRUE
niters = 0
while(notconverged){
	niters = niters + 1
	oldmu = Theta$mu
	
	W = E_step(as.matrix(D[,1:d]), Theta, Kmeans)
	
	Theta = M_step(as.matrix(D[,1:d]), W, Kmeans)
	print(Theta)
	
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
W = E_step(D[1:d], Theta, Kmeans, TRUE)
print(W)
#get true labels
clabels <- c("Iris-versicolor", "Iris-setosa", "Iris-virginica")
clabels = c(1,2,3)
cols = c("red", "blue", "green")
pchs = c(21, 22, 24)
TrueClusters <- list()
PredClusters <- list()

#get true and predicted (via EM) clusters
#pr = prcomp(D[,1:d])
#Z = pr$x[,1:2]
Z = D[,1:2]
for (i in 1:K){
	idx <- which(D[,d+1] == clabels[i])
	TrueClusters[[i]] <- idx
	idx <- which(Clusters == i)
	PredClusters[[i]] <- idx
}
print(PredClusters)

purity = compute_purity(TrueClusters, PredClusters)

plot(Z,pch=NA)
for (y in 1:3){
	C1 = which(D[,d+1] == clabels[y])
	for (i in 1:K){
		Zi = Z[C1[Clusters[C1]==i],]
		if (length(Zi) > 0){
			points(Zi, bg=cols[i], pch=pchs[y])
		}
	}
}
print(purity)

CK = matrix(0,4,K)
for (i in 1:K){
		Ci = which(Clusters==i)
		CK[1,i] = length(which(D[Ci,(d+1)]=="Iris-setosa"))
		CK[2,i] = length(which(D[Ci,(d+1)]=="Iris-virginica"))
		CK[3,i] = length(which(D[Ci,(d+1)]=="Iris-versicolor"))
		CK[4,i] = length(Ci)
}
print(CK)

#pr = prcomp(D)
#plot(pr$x[,1:2], pch=NA)
#cols = rainbow(K)
#for (i in 1:K){
	#Ci = which(Clusters==i)
	#points(pr$x[Ci,1:2], col=cols[i], pch=19)
#}


R = cbind(D,Clusters)
write.table(R, 'clusters.txt', row.names=FALSE, col.names=FALSE, sep=",")

#return(purity)
#}

##### MAIN #########
#datafile = 'data3c.txt'
#datafile = 'iris-PC.txt'
#K =3
#fullsigma=FALSE
#Kmeans=TRUE
#eps = 0.001 #stopping threshold

#D = read.table(datafile, sep=",")

#d = ncol(D)-1
#n = nrow(D)

#for (runs in 1:100){
	#SH <- sample(1:n)
	#Clusters <- c(rep(0,n))
	#pur = run_em()
	#print(c(runs,pur))
	#if (pur > 0.9){
		#print(SH)
		#break
	#}
#}
