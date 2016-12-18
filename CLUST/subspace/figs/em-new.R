#!/usr/bin/Rscript 
#setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/subspace/figs')
library(methods)
library(mvtnorm)
library(foreign)
library(MASS)

parse_args <- function(){
	Ieps = grep("eps", args)
	if (length(Ieps) > 0){
		eps <<- as.numeric(args[Ieps+1])
	}
	Iridge = grep("ridge", args)
	if (length(Iridge) > 0){
		ridgethresh <<- as.numeric(args[Iridge+1])
	}
	Iscale = grep("scale", args)
	if (length(Iscale) > 0){
		scaling <<- as.numeric(args[Iscale+1])
	}
	Idet = grep("det", args)
	if (length(Idet) > 0){
		detThreshold <<- as.numeric(args[Idet+1])
	}
	Ipr = grep("prcomp", args)
	if (length(Ipr) > 0){
		doprcomp <<- TRUE
		prcompDim <<- as.numeric(args[Ipr+1])
	}
	Isigma = grep("nofullsigma", args)
	if (length(Isigma) > 0){
		fullsigma <<- FALSE
	}
	Ifullafter = grep("fullafter", args)
	if (length(Ifullafter) > 0){
		fullafter <<- as.numeric(args[Ifullafter+1]) 
	}
	Ioracle = grep("oracle", args)
	if (length(Ioracle) > 0){
		oracle <<- TRUE
	}
	IKmeans = grep("^kmeans$", args)
	if (length(IKmeans) > 0){
		Kmeans <<- TRUE
	}
	Inoise = grep("noise", args)
	if (length(Inoise) > 0){
		noise <<- as.integer(args[Inoise+1])
	}
	Inorm = grep("normalize", args)
	if (length(Inorm) > 0){
		normalize <<- TRUE
		normalizeT <<- args[Inorm+1]
	}
	Ireldim = grep("reldim", args)
	if (length(Ireldim) > 0){
		relarg = args[Ireldim+1]
		ra = eval(parse(text=strsplit(relarg, ",")))
		reldim <<- unlist(lapply(ra, function(x){eval(parse(text=x))}))
		print(reldim)
	}
	Ipkmeans = grep("pkmeans", args)
	if (length(Ipkmeans) > 0){
		pkmeans <<- TRUE
	}
	IprojEM = grep("projEM", args)
	if (length(IprojEM) > 0){
		projEM <<- TRUE
	}
	Ialpha = grep("alpha", args)
	if (length(Ialpha) > 0){
		alpha <<- as.numeric(args[Ialpha+1])
	}
	Iprint = grep("print", args)
	if (length(Iprint) > 0){
		printparams <<- TRUE
	}
	Inoclass = grep("noclass", args)
	if (length(Inoclass) > 0){
		noclass <<- TRUE
	}
	Iplot = grep("plotfig", args)
	if (length(Iplot) > 0){
		plotfig <<- TRUE
	}
	Isep = grep("sep", args)
	if (length(Isep) > 0){
		sep <<- args[Isep+1]
	}
	Imaxiter = grep("maxiter", args)
	if (length(Imaxiter) > 0){
		maxiter <<- as.numeric(args[Imaxiter+1])
	}
}

choose_dim  <- function(vars){
	varD = sum(vars)
	cvar = cumsum(vars)
	mse = varD - cvar
	print(mse/mse[1])
	sa = which.min(mse/mse[1] <= alpha)
	m = mse[sa]/(sa-d)
	dists = c(rep(0,d))
	for (i in sa:d){
	    dists[i] = m*(i-d)-mse[i]
	}
	return (which.max(dists))
}


run_pkmeans <- function(D, Clusters, K, n, d){
	X <- D[,1:d]
	#last column will keep track of cluster assignment
	L <- Clusters
	X <- cbind(X,L) 

	not_converged <- TRUE
	SSEold <- 1
	iters <- 1
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
			dimvar <- svdi$d**2/ni
			dim_i  <- choose_dim(dimvar)
			svdi$v <- svdi$v[,1:dim_i]
			SVD_i[[i]] <- svdi

			#compute projection matrix
			pmi <- svdi$v %*% t(svdi$v)
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
				#minimize error
				XmM <- t(X[,1:d]) - M 
				P <- Q %*% XmM
				Sdist <- rowSums(t(P-XmM)*t(P-XmM))
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
		not_converged <- ifelse(abs(SSEnew-SSEold) > eps, TRUE, FALSE)
		cat(c("iters", iters, SSEnew, abs(SSEnew-SSEold), numChange, "\n"), sep=" ")
		iters <- iters+ 1
		SSEold <- SSEnew
	}
	
	L = X[,(d+1)]
	return (L)
}

apply_dmvnorm <- function(X, mu, S, k){
	EE = eigen(S, symmetric=TRUE, only.values=TRUE)
	detS = prod(EE$values)
	#print(c("detS", detS, k))
	if (detS < detThreshold){
		S = S + diag(ridgethresh,d,d)
		#print(c("newS", det(S)))
		res = dmvnorm(X,mu,S)
	}
	else{
		res = dmvnorm(X,mu,S)
	}
	return (res)
}


apply_projdist <- function(X, mu, S){
	EE = eigen(S, symmetric=TRUE)
	V = EE$vectors
	L = EE$values
	#dim_i  <- choose_dim(L)
	print(cumsum(L)/sum(L))
	print(alpha)
	dim_i <- which.min(cumsum(L)/sum(L) >= alpha)
	print(c(d-dim_i+1, dim_i))
	nV <- V[,1:dim_i]
	XmM <- t(X) - mu
	Px <- t(nV) %*% t(X)
	Px <- t(Px)
	print(dim(Px))
	newmu = c(rep(0,dim_i))
	if (dim_i == 1){
		newS = matrix(L[1:dim_i],1,1)
	} else {
		newS = diag(L[1:dim_i])
	}
	res <- dmvnorm(Px,newmu,newS)
	return (res)
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


###Expectation-Step
projE_step <- function(X, Theta, Assign=FALSE){
	W <- matrix(0, nrow=n, ncol=K)
	for (i in 1:K){
		m = Theta$mu[[i]]
		ss = Theta$sigma[[i]]
		if (d == 1){
			ss= as.matrix(Theta$sigma[[i]])
		}
		
		dmv = apply_projdist(X, m, ss)
		W[,i] = dmv*Theta$prior[i]
	}
	W = W/rowSums(W)
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
projM_step <- function(X, W){
	means = list()
	vars = list()
	priors <- colSums(W)/n #new prior
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
	
	return (list(mu=means, sigma=vars, prior=priors))
}



run_projEM <- function(D, Theta, K, n, d){
	X <- D[,1:d]
	#start proj EM iterations
	notconverged=TRUE
	niters = 0
	while(notconverged){
		niters = niters + 1
		oldmu = Theta$mu
		
		W = projE_step(X, Theta)
		
		Theta = projM_step(X, W)
		
		sumd = test_convergence (oldmu, Theta$mu)
		print(c("ITER", niters, sumd))
		if (sumd <= eps || niters > maxiter){
			notconverged = FALSE
		}
	}

	Clusters <- c(rep(0,n))
	W = projE_step(X, Theta, TRUE)
	return (Clusters)
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
datafile = args[1]
K = as.numeric(args[2])

eps = 0.0001 #stopping threshold
detThreshold = 0.000001
ridgethresh = 0.0001
scaling = 1
doprcomp = FALSE
fullsigma=TRUE
Kmeans=FALSE
maxiter = 500
noise = NULL
normalize = FALSE
normalizeT = "R" #range(R) or z(Z) normalization
oracle = FALSE
reldim = NULL
pkmeans = FALSE
projEM = FALSE
alpha = 0.7
printparams = FALSE
noclass = FALSE
plotfig = FALSE
fullafter = .Machine$integer.max
sep = ""

parse_args()

if (length(grep("arff", datafile)) > 0){
	X = read.arff(datafile)
} else {
	X = read.table(datafile,sep=sep)
}
n = nrow(X)
if (datafile == "test1.dat" || datafile == "test1.pca.dat"){
	TL = c(rep(0,n))
	np = 100
	for (i in 1:K){
		TL[((i-1)*np+1):(i*np)] <- c(rep(i,np))
	}
	X = cbind(X,TL) 
} else if (noclass){
	TL = c(rep(0,n)) #all points in same class
	X = cbind(X,TL)
}
d = ncol(X)-1
Z = X[,1:d]
if (normalize){
	if (normalizeT == "Z"){
		Z = scale(Z, center=TRUE, scale=TRUE)
	} else {
		Z = apply(Z, 2, function(x){(x-min(x))/(max(x)-min(x))})
	}
}
Z = Z*scaling
D = cbind(Z, X[,(d+1)])
print(c(n,d))

if (!is.null(reldim)){
	D = subset(D, select=c(reldim,(d+1)))
	d = length(reldim)
}

if (doprcomp){
	pr = prcomp(D[,1:d])
	D = cbind(pr$x[,1:prcompDim], D[,(d+1)])
	d = prcompDim
}

print(c("params", K, eps, scaling, ridgethresh, detThreshold, n, d))

#best centers found so far
use_fixed = FALSE
FM = list()
#FM[[1]] = c(0.379171203,  0.056546122, -0.003323449)
#FM[[2]] = c(0.53988846, 0.29295599, 0.03055884)
#FM[[3]] = c(0.34294427, 0.27967422, 0.05852964)

means = list()
sigmas = list()
priors = c(rep(1/K,K))
#initialize the clusters
if (oracle){
	ucl = unique(D[,(d+1)])
	K = length(ucl)
	for (i in 1:K){
		Dsub = D[D[,(d+1)] == ucl[i], ]
		ni = nrow(Dsub)
		priors[[i]] = ni
		means[[i]] = colMeans(Dsub[,1:d])
		sigmas[[i]] = cov(Dsub[,1:d])*(ni-1)/ni
	}
} else {
	notdone = TRUE
	siter = 1
	while(notdone){
		print(c("siter", siter, d))

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
}
#print(means)


Theta = list(mu=means,sigma=sigmas,prior=priors)
if (oracle){
	#tmpfullsigma = fullsigma
	#fullsigma = TRUE
	W = E_step(as.matrix(D[,1:d]), Theta, Kmeans)
	Theta = M_step(as.matrix(D[,1:d]), W, Kmeans)
	#fullsigma = tmpfullsigma
} else {
	W = E_step(as.matrix(D[,1:d]), Theta, TRUE)
	Theta = M_step(as.matrix(D[,1:d]), W, Kmeans)
}

#start EM iterations
notconverged=TRUE
niters = 0
while(notconverged){
	niters = niters + 1
	oldmu = Theta$mu
	
	if (niters >= fullafter){
		fullsigma = TRUE
	}

	W = E_step(as.matrix(D[,1:d]), Theta, Kmeans)
	
	Theta = M_step(as.matrix(D[,1:d]), W, Kmeans)
	
	sumd = test_convergence (oldmu, Theta$mu)
	print(c("ITER", niters, sumd))
	if (sumd <= eps || niters > maxiter){
		notconverged = FALSE
	}
}

if (printparams){
	print ('Cluster Parameters')
	#print(Theta)
	print(Theta$prior)
	print(Theta$mu)
	if (!Kmeans){
		for (i in 1:K){
			if (fullsigma){
				print(Theta$sigma[[i]])
			} else {
				print(diag(Theta$sigma[[i]]))
			}
		}
	}
}
print(c('Number of Iterations', niters))

#get final clusters
Clusters <- c(rep(0,n))
W = E_step(D[,1:d], Theta, Kmeans, TRUE)

if (pkmeans){
	Clusters = run_pkmeans(D,Clusters,K,n,d)
} 


#check how many match
CT <- table(Clusters, D[,(d+1)])

print(CT)
if (!is.null(noise)){
	nfval <- 0.0
	elimrow  <- NULL
	mn = which.max(CT[,noise])
	if (CT[mn,noise]*1.0/sum(CT[mn,]) >= 0.5){
		elimrow <- mn
	}
	print(c("noise", noise, elimrow))
	for (i in 1:nrow(CT)){
		if (!is.null(elimrow) && i == elimrow){ next}
		mC = max(CT[i,-noise])
		mj = which.max(CT[i,-noise])
		prec  <- mC / sum(CT[i,])
		recall <- mC / sum(CT[,mj])
		fvali  <- (2.0*prec*recall)/(prec+recall)
		nfval  <- nfval + fvali
	}
	nfval = nfval/(K-1)
	cat('nF', nfval, '\n')
}

purity <- 0.0
fval <- 0.0
entropy <- 0.0
for (i in 1:nrow(CT)){
	ei <- 0.0
	for (j in 1:ncol(CT)){
	   pij = CT[i,j]*1.0/sum(CT[i,])
	   if (pij > 0){
		ei = ei - (pij * log2(pij))
	   }
	}
	entropy = entropy + sum(CT[i,])*ei/n
	mC = max(CT[i,])
	mj = which.max(CT[i,])
    purity <- purity + mC
	prec  <- mC / sum(CT[i,])
	recall <- mC / sum(CT[,mj])
	fvali  <- (2.0*prec*recall)/(prec+recall)
	#print(c("i", mC, mj, prec, recall, sum(CT[i,]), sum(CT[,mj]), fvali))
	fval  <- fval + fvali
}
fval = fval/K

#diagonal purity, max matching
dPurity = 0.0
nCT = ncol(CT)
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
cat('F', fval, '\n')
cat('Purity', purity, '\n')
cat ('dPurity', dPurity, '\n')
cat ('entropy', entropy, '\n')

if (plotfig){
	colary = sample(colours(), K)
	plot(X[,1],X[,2])
	for (i in 1:K){
		Pi = which(Clusters == i)
		points(X[Pi,1], X[Pi,2], pch=19, col=colary[i])
	}
}

if (projEM){
	Clusters = run_projEM(D,Theta,K,n,d)
}
#check how many match
#CT <- table(D[,(d+1)],Clusters)
#print(CT)
