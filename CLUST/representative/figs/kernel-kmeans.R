setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/representative/figs')
library(igraph)
library(Cairo)

phi2h <- function(x){
	r = c(rep(0,((d+1)*(d+2)/2)-1));
	k=1;
	for (i in 1:d){
		r[k] = sqrt(2)*x[i]
		k = k+1
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


run_kmeans <- function(D,K,d,n){
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
	niters = 0
	notconverged = TRUE
	while(notconverged){
		niters = niters + 1
		oldmu = means
	
		#E_step
		W <- matrix(0, nrow=n, ncol=K)
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
	
	W <- matrix(0, nrow=n, ncol=K)
	M = matrix(unlist(means),d,K)
	dd = matrix(0,n,K)
	for (i in 1:K){
		dd[,i] = colSums((t(D)-M[,i])^2)
	}
	for (j in 1:n){
		ci = which.min(dd[j,])
		W[j,ci] = 1
	}
	return (W)
}


fnames = c("iris.txt", "data3c.txt", "kerneldata.txt")
fname = fnames[3]

X = read.table(fname, sep=",")
d = ncol(X)-1
D = as.matrix(X[,1:d])
K = 3

n = nrow(D)
M = matrix(0,n,n)

runkmeans=0

kernelvals = c("gaussian", "gaussian-local", "linear", "quadratic")

kernel = kernelvals[1]

if (kernel == "gaussian"){ #Gaussian
	for (i in 1:n){
		Zi = t(t(D) - D[i,]);
		rsi = rowSums(Zi*Zi);
		sigma=1.5
		di = exp(-1*rsi/(2*sigma^2));
		M[i,] = di;
	}
} else if (kernel == "gaussian-local"){ #Gaussian Local
	#from Self-Tuning Spectral Clustering; Zelnik-Manor/Perona
	delmat = matrix(0,n,n)
	sigma = c(rep(0,n))
	knn=10
	for (i in 1:n){
		Zi = t(t(D) - D[i,]);
		delmat[i,] = rowSums(Zi*Zi);
		Oi = order(delmat[i,])
		sigma[i] = delmat[i,Oi[knn]]
	}
	for (i in 1:n){
		di = exp(-1*delmat[i,]/(sqrt(sigma*sigma[i])));
		M[i,] = di;
	}
} else if (kernel == "linear"){ #Linear
	bias = 0
	M = D %*% t(D) + bias
} else if (kernel == "quadratic"){ #Quadratic
	bias = 0
	M = (D %*% t(D) + bias)^2
} 

#m = mean(M)
#s = sqrt(var(c(M)))
#M[M < m+1*s] = 0
#M[M < 0.000001] = 0


#parition points into k clusters randomly
bs = ceiling(n/K)
SH <- sample(1:n)
W = matrix(0,n,0)
for(i in 1:K){
        Pi = c(rep(0, n))
        l = (i-1)*bs+1
		u = min(n,i*bs)
        Pi[SH[l:u]] = 1.0
		W = cbind(W,Pi)
}

if (runkmeans){
	W = run_kmeans(D,K,d,n)
}
#print(W)
#print(colSums(W))

#start Kmeans iterations
notconverged=TRUE
niters = as.integer(0)

oldSSE = 0
eps = 0.001
oneK = matrix(1,K,1)
oneN = matrix(1,n,1)
while(notconverged){
	Z = matrix(0,n,K)
	SSE = 0
	for (i in 1:K){
		nz = which(W[,i] > 0)
		ni = length(nz)
		R = M[nz,nz]
		cK = sum(R)/ni
		SSE = SSE + sum(diag(R)) - cK
		S = (2.0/ni)*colSums(M[nz,])
		Z[,i] = cK/ni-S
	}
	
	#Wbn = t(t(W)/colSums(W))
	#S = 2*(M %*% Wbn)
	#T = oneN %*% diag(t(Wbn) %*% M %*% Wbn)
	#Z2 =  T-S
	
	AssignVec = apply(Z,1,which.min)
	nW = matrix(0,n,0)
	nchanges = 0
	for (i in 1:K){
		Pi = c(rep(0, n))
		Pi[which(AssignVec == i)] = 1.0
		nW = cbind(nW,Pi)
		isect = intersect(which(W[,i] == 1),which(nW[,i]==1))
		nchanges = nchanges + length(isect)
	}
	diff = (n-nchanges)/n
	#diff = abs(SSE-oldSSE)
	if (diff <= eps || niters > 50) {notconverged = FALSE}
	W = nW
	oldSSE = SSE
	print(c("ITER", niters, diff, SSE, nchanges))
	niters = niters + 1
}

colSums(W)

#compute SSE value
SSE = 0
for (i in 1:K){
	nz = which(W[,i] > 0)
	ni = length(nz)
	SSE = SSE + (sum(diag(M)[nz]) - sum(M[nz,nz])/ni)
}
SSE

if (fname == fnames[1]){
	pr = prcomp(D)
	plot(pr$x[,1:2], pch=NA)
	clas = c('Iris-setosa','Iris-virginica','Iris-versicolor')
	cols = c('green', 'red', 'blue')
	pchs = c(22, 24, 21)

	for (y in 1:3){
		C1 = which(X[,5] == clas[y])
		for (i in 1:K){
			Di = pr$x[C1[AssignVec[C1] == i],1:2]
			if (length(Di) > 0){
				points(Di, bg=cols[i], pch=pchs[y])
			}
		}
	}

	CK = matrix(0,4,K)
	for (i in 1:K){
			Ci = which(AssignVec==i)
			CK[1,i] = length(which(X[Ci,5]=="Iris-setosa"))
			CK[2,i] = length(which(X[Ci,5]=="Iris-virginica"))
			CK[3,i] = length(which(X[Ci,5]=="Iris-versicolor"))
			CK[4,i] = length(Ci)
	}
	print(CK)
} else {
	plot(D, pch=NA)
	cols = rainbow(K)
	for (i in 1:K){
		Ci = which(AssignVec==i)
		points(D[Ci,], col=cols[i], pch=19)
	}

}
print(c("SSE", SSE))
R = cbind(X,AssignVec)
write.table(R, 'nlclusters.txt', row.names=FALSE, col.names=FALSE, sep=",")

#if (kernel == "quadratic"){
	#phiX = t(apply(D, 1, phi2h))
	#pr = prcomp(phiX)
	#plot(pr$x[,1:2])
#}

#get true labels
clabels = c(1,2,3)
#clabels <- c("Iris-versicolor", "Iris-setosa", "Iris-virginica")
cols = c("red", "blue", "green")
pchs = c(21, 22, 24)
TrueClusters <- list()
PredClusters <- list()

#get true and predicted (via EM) clusters
#pr = prcomp(D[,1:d])
#Z = pr$x[,1:2]
for (i in 1:K){
	idx <- which(X[,d+1] == clabels[i])
	TrueClusters[[i]] <- idx
	idx <- which(AssignVec == i)
	PredClusters[[i]] <- idx
}

purity = compute_purity(TrueClusters, PredClusters)
print(purity)
