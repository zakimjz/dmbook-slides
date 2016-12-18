setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/spectral/figs')
library(igraph)
library(Cairo)

compute_purity <-function(T, P){
	summ = 0
	Tc = max(T) #num true classes
	Tp = max(P) #num pred classes
	for (i in 1:Tp){
		Cp = which(P == i)
		isect = c(rep(0,Tc))
		for (j in 1:Tc){
			Ct = which(T == j)
			isect[j] = length(intersect(Ct,Cp))
		}
		summ = summ + max(isect)
	}
	return(summ/n)
}


run_kmeans <- function(D,K,d,n,attrW=c()){
	eps = 0.001
	means = list()
	#attrW = seq(d,1) #weighted K-means
	#attrW = attrW/sum(attrW)
	for (i in 1:K){
		mu = c(rep(0,d))
		for (j in 1:d){
			r  =range(D[,j])
			mu[j] = runif(1,r[1], r[2])
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
			if (length(attrW) > 0){
				dd[,i] = colSums(((t(D)-M[,i])^2)*attrW)
			} else {
				dd[,i] = colSums((t(D)-M[,i])^2)
			}
		}
		for (j in 1:n){
			ci = which.min(dd[j,])
			W[j,ci] = 1
		}
		
		#M_step
		means = list()
		for (i in 1:K){
			if (sum(W[,i]) > 0){
				mi = colSums(D*W[,i])/sum(W[,i])
			} else {
				mi = c(rep(0,d))
				for (j in 1:d){
					r  =range(D[,j])
					mi[j] = runif(1,r[1], r[2])
				}
			}
			means[[i]] <- mi
		}

		sumd = 0
		notconverged = TRUE
		for (i in 1:K){
				di = oldmu[[i]]-means[[i]]
				diff = di * di #magnitude of diff
				sumd <- sumd + sum(diff)
		}

		#browser()
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
	clust = c(rep(0,n))
	for (j in 1:n){
		ci = which.min(dd[j,])
		W[j,ci] = 1
		clust[j] = ci
	}
	print(means)
	return (list(cluster=clust))
}

run_mcl  <- function(A,p){
	cs = rowSums(A)
	Z = A/cs
	oZ = Z;
	eps = 0.001
	iter = as.integer(1)
	diff = 1
	while(diff > eps){
		Z = Z%*%Z
		Z = Z^p
		cs = rowSums(Z)
		Z = Z/cs
		diff = sum((Z-oZ)^2)
		print(c(iter, diff))
		oZ = Z
		iter = iter+1
	}

	Z[which(Z<0.001)] = 0 #make very small values 0
	Z[which(Z>0)]=1

	zg = graph.adjacency(Z,mode="undirected")
	clust = clusters(zg)

	return (list(cluster=(clust$membership+1),Z=Z))
}


####MAIN####
graphconstalgs = c("full", "knn", "mutualknn")
distfuncs = c("reciprocal", "exp")
clusteralgs = c("normalized", "ratio", "modularity", 
		"average", "avgnorm", "MCL")

alg=graphconstalgs[3]
binarymatrix = FALSE #for use with full
distfunc = distfuncs[2]
clusteralg = clusteralgs[1]
symmetric=FALSE #for normalized laplacian

print(c(alg, distfunc, clusteralg, symmetric))
knn=16
k = 3 #num clusters
p = 2.5 #powers of the elements for MCL

graph="iris"
#graph="example"
#graph = "karate"

if (graph == "iris"){
	X = read.table('iris.txt', sep=",")
	D = as.matrix(X[1:4])
	n = nrow(D);
	M = matrix(0,n,n);

	for (i in 1:n){
		Zi = t(t(D) - D[i,]);
		rsi = rowSums(Zi*Zi);
		if (distfunc == "exp"){
			sigma=1
			rsi = rsi/(2*sigma^2)
			di = exp(-1*rsi);
		} else if (distfunc == "reciprocal"){
			di = sqrt(rsi);
		}
		M[i,] = di;
	}
	if (distfunc == "reciprocal"){
		maxv = max(M)
		minv = min(M)
		M = maxv-M
		M = M/(maxv-minv)
	}


	g = graph.empty(n=n, directed=FALSE)
	CC = c(rep(0,n))

	for (i in 1:n){
		if (X[i,5] == "Iris-setosa") {CC[i] = "red";}
		if (X[i,5] == "Iris-versicolor") {CC[i] = "blue";}
		if (X[i,5] == "Iris-virginica") {CC[i] = "green";}
	}
	V(g)$color = CC
	V(g)$names = c(as.integer(X[,5]))
	V(g)$clust = c(rep(1,n))

	if (alg == "full"){
		Z = M[M>0]
		m = mean(Z)
		s= sd(Z)
		simthresh = m + 1.75*s
		for (i in 1:n){
			for (j in i:n){
				if (M[i,j] > simthresh){
					if (binarymatrix == TRUE){
						#g = add.edges(g, c(i-1,j-1), weight=1.0)
						g[i,j,attr="weight"] = 1.0
					} else {
						#g = add.edges(g, c(i-1,j-1), weight=M[i,j])
						g[i,j,attr="weight"] = M[i,j] 
					}
				}
			}
		}
	} else if (alg == "knn"){
		MKNN = list()
		for (i in 1:n){
			Oi = order(M[i,], decreasing=TRUE)
			for (j in knn:n){
				if (M[i,Oi[j]] != M[i,Oi[knn]]) {
					break
				}
			}
			MKNN[[i]] = Oi[1:(j-1)]
		}

		for (i in 1:n){
			for (j in MKNN[[i]]){
				#if (i != j){
					#g = add.edges(g, c(i-1,j-1), weight=M[i,j])
					g[i,j,attr="weight"] = M[i,j]
					#}
			}
		}
	} else if (alg == "mutualknn"){
		MKNN = list()
		for (i in 1:n){
			Oi = order(M[i,], decreasing=TRUE)
			for (j in knn:n){
				if (M[i,Oi[j]] != M[i,Oi[knn]]) {
					break
				}
			}
			MKNN[[i]] = Oi[1:(j-1)]
		}
		for (i in 1:n){
			for (j in MKNN[[i]]){
				if (i %in% MKNN[[j]]){# && i != j){
					#g = add.edges(g, c((i-1),(j-1)), weight=M[i,j])
					g[i,j,attr="weight"] = M[i,j]
				}
			}
		}
	}
} else if (graph == "example"){
	n=7
	M = matrix(c(0,1,0,1,0,1,0, 
				 1,0,1,1,0,0,0, 
				 0,1,0,1,0,0,1, 
				 1,1,1,0,1,0,0, 
				 0,0,0,1,0,1,1, 
				 1,0,0,0,1,0,1, 
				 0,0,1,0,1,1,0), n,n)
	I = diag(n)
	if (clusteralg == "MCL"){
		M = M+I
	}

	g = graph.empty(n=n, directed=FALSE)
	V(g)$color=c(rep("red",4), rep("blue",3))
	V(g)$names=c(1,1,1,1,2,2,2)
	V(g)$clust = c(rep(1,n))


	#n=4
	#M = matrix(c(2.33, 1.65, 1.87, 2.18, 
	#			 1.65, 1.62, 1.69, 1.58,
	#			 1.87, 1.69, 1.81, 1.78,
	#			 2.18, 1.58, 1.78, 2.05
	#			 ), n,n)
	#M = matrix(c(0,1,0,1, 
	#			 1,0,1,1, 
	#			 0,1,0,1, 
	#			 1,1,1,0), n,n)


	#g = graph.empty(n=n, directed=FALSE)
	#V(g)$color=c(rep("red",4))
	#V(g)$names=c(1,1,1,1)
	#V(g)$clust = c(rep(1,n))


	for (i in 1:(n-1)){
		for (j in (i+1):n){
			if (M[i,j] > 0){
				#g = add.edges(g,c(i-1,j-1),weight=M[i,j])
				g[i,j,attr="weight"] = M[i,j]
			}
		}
	}
	print(M)
	print(g)
} else if (graph == "karate"){
	M = as.matrix(read.table("karate.txt"))
	n = nrow(M)
	I = diag(n)
	if (clusteralg == "MCL"){
		M = M+I
	}
	g = graph.empty(n=n, directed=FALSE)
	V(g)$clust = c(rep(1,n))
	V(g)$names = c(rep(1,n))
	V(g)$color=c(rep("red",n))
	for (i in 1:(n-1)){
		for (j in (i+1):n){
			if (M[i,j] > 0){
				#g = add.edges(g,c(i-1,j-1),weight=M[i,j])
				g[i,j,attr="weight"] = M[i,j]
			}
		}
	}
	print(g)
}
#add all missing links between components
comps = clusters(g)
print(comps)
if (comps$no > 1){
	for (i in 1:(comps$no-1)){
		Vi = which(comps$membership == i)
		for (j in (i+1):comps$no){
			Vj = which(comps$membership == j)
			OO = sort(M[Vi,Vj], decreasing=TRUE)
			knnv = OO[min(length(OO),knn)]
			for (a in Vi){
				for (b in Vj){
					if (M[a,b] >= knnv){
						#print(c(i,j,a-1,b-1,M[a,b], knnv))
						if (alg != "full" || !binarymatrix){
							#g = add.edges(g, c(a-1,b-1), weight=M[a,b])
							g[a,b,attr="weight"] = M[a,b]
						}
					}
				}
			}
		}
	}
}

E(g)$deleted = c(rep(0,ecount(g)))

#delete all small components, and keep only the largest one
comps = clusters(g)
print(comps)
#del_c = order(comps$csize, decreasing=TRUE)[3:comps$no]
#g=delete.vertices(g,which(comps$membership %in% (del_c-1))-1)

if (graph == "iris" && file.exists('iris-base-layout.txt')){
	l = as.matrix(read.table('iris-base-layout.txt'))
} else {
l = layout.fruchterman.reingold(g, repulserad=225, niters=2000,
		area=500)
l = layout.norm(l,0,5,0,5,0)
}

if (graph == "iris"){
	A = matrix(0,n,n)
	#print(ecount(g))
	#print(get.edges(g,1:ecount(g)))
	for (i in 1:ecount(g)){
		e = get.edge(g,i)
		A[e[1],e[2]] = E(g)$weight[i]
		A[e[2],e[1]] = E(g)$weight[i]
	}
} else if (graph == "example" || graph == "karate"){
	A = M
}

plot(g,layout=l, vertex.size=5, vertex.label=NA)

if (clusteralg != "MCL"){
	#spectral clustering
	D = diag(rowSums(A))
	if (clusteralg == "ratio"){
		L = D-A
	} else if (clusteralg == "normalized") {
		if (symmetric){
			ID = diag(diag(D^{-1/2}))
			L = ID %*% (D-A) %*% ID
		} else {
			ID = diag(diag(D^{-1}))
			L = ID%*% (D-A)
		}
	} else if (clusteralg == "modularity"){
		d = as.matrix(rowSums(A),n,1)
		sd = sum(d)
		E = d%*%t(d)
		E = E/sd
		L = (A-E)/sd
	} else if (clusteralg == "average"){
		L = A
	} else if (clusteralg == "avgnorm"){
		ID = diag(diag(D^{-1}))
		L = ID %*% A
	}
	#print(D)
	#print (L)
	ss=eigen(L)
	print(ss$values[(n):(n-k)])
	print(ss$vectors[,(n):(n-k)])
	n = vcount(g)

	if (clusteralg == "normalized" || 
		clusteralg == "ratio"){
		#get the smallest eigenvectors
		U = ss$vectors[,(n):(n-k+1)]
	} else if (clusteralg == "modularity" || 
				clusteralg == "average"){
		#get the largest eigenvectors
		U = ss$vectors[,1:k]
	} else if (clusteralg == "avgnorm"){
		Ov = order(ss$values, decreasing=TRUE)
		U = ss$vectors[,Ov[1:k]]
	}


	attrW = c()
	if (clusteralg == "modularity" || clusteralg == "avgnorm"){	
		attrW = seq(k,1) #weighted K-means
		attrW = attrW/sum(attrW)
		#attrW = ss$values[1:k]/sum(ss$values[1:k])
	}
  
	#if (clusteralg == "normalized" && symmetric){
	for (i in 1:n){
		if (sum(U[i,]*U[i,]) > 0){
		  U[i,] = U[i,]/sqrt(sum(U[i,]*U[i,]))
		}
	}
	#}
	print(U)
	write.table(U, 'Y.txt', row.names=FALSE, col.names=FALSE)
	#kk = kmeans(U,k)
	
	kk = run_kmeans(U,k,k,n,attrW)
} else if (clusteralg == "MCL"){
	kk = run_mcl(A,p)
}
print(kk)
k = max(kk$cluster) #re-adjust clusters if fewer/more are returned

for (i in 1:k){
	V(g)$clust[which(kk$cluster == i)]=i
}
print(g)
#compute normalized cut value
Jw = 0
for (i in 1:k){
	Ci = which(kk$cluster==i)
	voli = sum(A[Ci,])
	cuti = sum(A[Ci,-Ci])
	Jw = Jw + (cuti/voli)
}
print(c("normalized cut value", Jw))

#compute the modularity value
Jmod=0
d = as.matrix(rowSums(A),n,1)
sd = sum(d)
for (i in 1:k){
	Ci = which(kk$cluster==i)
	Ai = A[Ci,Ci]
	Di = d[Ci]%*%t(d[Ci])/sd
	modi = sum(Ai-Di)
	Jmod = Jmod+modi
}
Jmod = Jmod/sd

CutVal = 0
ng = g #new graph

for (i in 1:ecount(g)){
	e = get.edge(g,(i))
	u = e[1]
	v = e[2]
	if (V(g)$clust[u] != V(g)$clust[v]){
		wt = E(g)$weight[i]
		CutVal = CutVal + wt
		ng = delete.edges(ng,E(ng,P=c(u,v)))
		E(g)$deleted[i]=1
	}
}
print(c("cut value", CutVal))

colary = c("red","green","blue")
V(ng)$color = colary[kk$cluster]


plot(ng,layout=l, vertex.size=5, vertex.label=NA)

labs = c("Iris-setosa", "Iris-versicolor", "Iris-virginica")
if (graph == "iris"){
	cat (file='clusters.txt')

	if (clusteralg == "MCL"){
		for (v in V(g)){
			att=0 #is the vertex an attractor
			if (kk$Z[v,v] > 0){
				att=1
			}
			cat (file='clusters.txt', append=TRUE, 
				"v", v, l[v,1],
				l[v,2], labs[V(g)$names[v]], 
				V(g)$clust[v], att, "\n")
		}
		#write the graph edges	
		for (i in 1:ecount(g)){
			e = get.edge(g,(i))
			cat (file='clusters.txt',  append=TRUE, 
				"e", e[1], e[2], E(g)$weight[i], E(g)$deleted[i], "\n")
		}
		#write the attractor edges
		for (i in 1:n){
			for (j in 1:n){
				if (kk$Z[i,j] > 0 && i!=j){
					cat (file='clusters.txt',  append=TRUE, 
					"de", i, j, "\n")
				}
			}
		}
	} else{
		for (v in V(g)){
			cat (file='clusters.txt', append=TRUE, 
				"v", v, l[v,1],
				l[v,2], labs[V(g)$names[v]], 
				V(g)$clust[v], "\n")
		}
		for (i in 1:ecount(g)){
			e = get.edge(g,i)
			cat (file='clusters.txt',  append=TRUE, 
				"e", e[1], e[2], E(g)$weight[i], E(g)$deleted[i], "\n")
		}
	}
}

if (graph == "iris"){
	CK = matrix(0,4,k)
	for (i in 1:k){
		Ci = which(kk$cluster==i)
		CK[1,i] = length(which(labs[V(g)$names[Ci]]=="Iris-setosa"))
		CK[2,i] = length(which(labs[V(g)$names[Ci]]=="Iris-virginica"))
		CK[3,i] = length(which(labs[V(g)$names[Ci]]=="Iris-versicolor"))
		CK[4,i] = length(Ci)
	}
	print(CK)
}

purity = compute_purity(V(g)$names, kk$cluster)
print(c("cut value, Jw, purity, Jmod", CutVal, Jw, purity, Jmod))
