setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/subspace/figs')

# definition of a function for parameterized data simulation
sim.orclus <- function(k = 3, nk = 1, d = 1, l = 4, 
					   sd.cl = 0.5, sd.rest = 1, locshift = 1){
	### input parameters for data generation
	# k	number of clusters
	# nk observations per cluster
	# d original dimension of the data
	# l subspace dimension where the clusters are concentrated
	# sd.cl (univariate) standard deviations for data generation (within cluster subspace concentration)
	# sd.rest standard deviations in the remaining space
	# locshift parameter of a uniform distribution to sample different cluster means

	x <- NULL
	for(i in 1:k){
		# cluster centers
		apts <- locshift*matrix(runif(l*k), ncol = l)
		# sample points in original space
		xi.original <- cbind(matrix(rnorm(nk * l, sd = sd.cl), ncol=l) 
							 + 
							 matrix(rep(apts[i,], nk), ncol = l, 
									byrow = TRUE),

		matrix(rnorm(nk * (d-l), sd = sd.rest), ncol = (d-l)))
		# subspace generation
		sym.mat <- matrix(nrow=d, ncol=d)
		for(m in 1:d){
			for(n in 1:m){
			sym.mat[m,n] <- sym.mat[n,m] <- runif(1)
			}
		}
		subspace <- eigen(sym.mat)$vectors
		# transformation
		xi.transformed <- xi.original %*% subspace
		x <- rbind(x, xi.transformed)
	}
	clids <- rep(1:k, each = nk)
	result <- list(x = x, cluster = clids)
	return(result)
}

D <- sim.orclus(k = 3, nk = 100, d = 3, l = 2, 
					  sd.cl = 0.05, sd.rest = 1, locshift = 1)

W = cbind(D$x, D$cluster)
write.table(W, file="simdata1.txt", row.names=FALSE, col.names=FALSE)
library(orclus)
res = orclus(x=W[,1:3], k=3, l=3, k0=20, a=0.75)
CT <- table(res$cluster, W[,4])
purity <- 0
for (i in 1:ncol(CT)){
 	purity <- purity + max(CT[,i])
}
purity <- purity/300

