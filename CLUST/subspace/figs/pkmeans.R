#!/usr/bin/env Rscript
require(stats)

kmeans <- function(D, means, K,n,d){
	eps = 0.001
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

initial_partition  <- function(D,K,n,d){
	if (select_means) {
		notdone = TRUE
		siter = 1
		while(notdone){
			print(c("siter", siter))
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

			M = matrix(unlist(means),d,K)
			dd = matrix(0,n,K)
			for (i in 1:K){
				dd[,i] = colSums((t(D)-M[,i])^2)
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
		
		if (run_kmeans){
			L <- kmeans(D, means, K, n, d)
		}
	} else {
		L <- sample(1:3, n, replace=TRUE) #random point assignment
		print(c("random seed", .Random.seed))
	}
	return (L)
}


choose_dim  <- function(vars){
	varD = sum(vars)
	cvar = cumsum(vars)
	mse = varD - cvar
	sa = min(which(mse/mse[1] <= alpha))
	m = mse[sa]/(sa-d)
	dists = c(rep(0,d))
	for (i in sa:d){
	    dists[i] = m*(i-d)-mse[i]
	}
	return (which.max(dists))
}

####MAIN###
args <- commandArgs(TRUE)
setwd("/Users/zaki/research/DataMiningBook/dm08/CLUST/subspace/figs")
#datafile = "test1.dat"
#datafile = "test1.pca.dat"
#datafile = "iris-PC2.txt"
datafile = args[1]
K <- as.numeric(args[2])
alpha <- 0.2
var_threshold = 0.2
if (length(args) > 2){
  var_threshold = as.numeric(args[3])
}
alpha = var_threshold

scale = 1
D <- read.table(datafile)*scale
d <- ncol(D)
n <- nrow(D)
if (datafile == "test1.dat" || datafile == "test1.pca.dat"){
	TL = c()
	for (i in 1:K){
		TL <- c(TL, rep(i,100))
	}
	D <- cbind(D, TL) #attach true cluster labels
} else{
	d = d-1 #last dim is class
}
X <- D[,1:d]
changeThresh = 10
tol <- 0.0001
varmax =  TRUE
run_kmeans = FALSE
select_means = TRUE

#best centers found so far
use_fixed = FALSE
FM = list()
FM[[1]] = c(0.3301256, 1.1180706, 4.6055505)
FM[[2]] = c(7.8478957, 0.7246735, 4.2788653)
FM[[3]] = c(1.567903, 1.440511, 5.204630)

#Initial partition
L <- initial_partition(X,K,n,d)

#last column will keep track of cluster assignment
X <- cbind(X,L) 

not_converged <- TRUE
SSEold <- 1
iters <- 1
stuckcnt  <- 0
numChange  <- 0
while(not_converged){
	SVD_i <- c() #keep svd dimension info
	Mean_i <- c() #shift for correct distance
	Projection_i <- c() #projection matrices
	#Compute new subspaces
	Valid_i  <- c(rep(TRUE,K))
	prevL <- L
	prevnumChange <- numChange
	for (i in 1:K){
		Di <- subset(X, X[(d+1)] == i)
		ni <- nrow(Di)
		if (ni == 0) {
			print(c("not VALID", i))
			Valid_i[i] = FALSE
			next
		}
		Di <- Di[,1:d]
		Di = as.matrix(Di)
		Mi <- colSums(Di)/ni
		DMi <- t(t(Di)-Mi)
		svdi <- svd(DMi, nu=0) #not interested in U
		#find subspace dimensionality
		if (varmax){
			#print(svdi$d)
			dimvar <- svdi$d**2/ni
			#dimrat <- cumsum(dimvar)/sum(dimvar)
			#dim_i <- min(which(dimrat >= var_threshold))
			dim_i  <- choose_dim(dimvar)
			#select relevant dimensions
			svdi$v <- svdi$v[,1:dim_i]
		} else{
			#print(svdi$d)
			dimvar <- rev(svdi$d**2/ni)
			dimrat <- cumsum(dimvar)/sum(dimvar)

			wcc = which(dimrat <= var_threshold)
			if (length(wcc) > 0){
				midx = max(wcc)
			} else {
				midx = d #all dims are important
			}
			dim_i <- midx
			#dim_i <- max(which(dimrat <= var_threshold))
			#print(c("dimensionality", i, dim_i, dimrat))
			#select relevant dimensions
			svdi$v <- svdi$v[,(d-dim_i+1):d]
		}
		SVD_i[[i]] <- svdi

		#compute projection matrix
		pmi <- svdi$v %*% t(svdi$v)
		#print(pmi)
		Projection_i[[i]] <- pmi
		Mean_i[[i]] <- Mi
	}

	#Assign point to subspaces, to create new partition
	L = c(rep(0,n)) #new cluster assignment vector
	Distances <- matrix(1000000,n,K) #compute distance of each point to subspace
	for (i in 1:K){
		if (Valid_i[i]){
			M = Mean_i[[i]]
			Q = Projection_i[[i]]
			if (varmax){
				#minimize error
				XmM <- t(X[,1:d]) - M 
				P <- Q %*% XmM
				Sdist <- rowSums(t(P-XmM)*t(P-XmM))
			} else{
				#minimize SSE
				pM =  Q %*% M
				dim(pM) = NULL #make it into a regular list
				Px = Q %*% t(X[,1:d])
				Sdist <- rowSums(t(Px-pM)**2)
			}
			Distances[,i] <- Sdist
		}
	}	
	L <- apply(Distances, 1, which.min)
	MinDist <- apply(Distances, 1, min)
	X[,(d+1)] <- L #replace assignment column with new one
	CT <- table(L,prevL)
	numChange = n-sum(diag(CT))
	#test convergence
	SSEnew <- sum(MinDist)
	not_converged <- ifelse(abs(SSEnew-SSEold) > tol, TRUE, FALSE)
	#not_converged  <- ifelse(numChange <= changeThresh, FALSE, TRUE)
	#if (numChange == prevnumChange){ 
	#	stuckcnt = stuckcnt + 1
	#}
	#if (not_converged){
	#	not_converged  <- ifelse(stuckcnt > 10, FALSE, TRUE)
	#}
	cat(c("iters", iters, SSEnew, abs(SSEnew-SSEold), numChange, "\n"), sep=" ")
	iters <- iters+ 1
	SSEold <- SSEnew
}

#print cluster information
for (i in 1:K){
	Di <- as.matrix(subset(X, X[(d+1)] == i))
	Di <- Di[,1:d]
	ni  <- nrow(Di)
	if (ni == 0) {next}
	print(c("cluster", i, ": size=", ni))
	Mi <- colSums(Di)/ni
	DMi <- t(t(Di)-Mi)
	svdi <- svd(DMi, nu=0)
	vars = svdi$d**2/ni
	if (varmax){
		#cc <- cumsum(vars)/sum(vars)
		#dim_i = min(which(cc >= var_threshold))
		dim_i = choose_dim (vars)
	} else {
		vars = rev(vars) #reverse the vals
		cc <- cumsum(vars)/sum(vars)
		wcc = which(cc <= var_threshold)
		if (length(wcc) > 0){
			midx = max(wcc)
		} else {
			midx = d #all dims are important
		}
		dim_i <- midx
	}
	print(vars)
	#print(cc)
	print(c("dimensionality", dim_i))
	print(Mi)
	print(cov(Di)*(ni-1)/ni)
	print(svdi$v)
}


#check how many match
TP <- cbind(D[,(d+1)], X[,(d+1)])
CT <- table(TP[,1],TP[,2])

print(CT)
purity <- 0
for (i in 1:ncol(CT)){
	purity <- purity + max(CT[,i])
}
print(c(purity,n))
purity <- purity/n
cat('Purity', purity, '\n')
