setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/spectral/figs')
library(igraph)
library(Cairo)
alg = "full"
#alg = "knn"
knn=5
X = read.table('iris.txt', sep=",")
D = as.matrix(X[1:4])
M = matrix(0,150,150);

for (i in 1:150){
	Zi = t(t(D) - D[i,]);
	rsi = rowSums(Zi*Zi);
	di = exp(-1*rsi);
	M[i,] = di;
}

Z = M[M>0]
m = mean(Z)
s= sd(Z)

simthresh = m + 1.75*s
simthresh2 = m + 1*s
simthreshKNN = m + 2*s

g = graph.empty(n=151, directed=FALSE)
CC = c(rep(0,151))

for (i in 1:150){
	if (X[i,5] == "Iris-setosa") {CC[i+1] = "red";}
	if (X[i,5] == "Iris-versicolor") {CC[i+1] = "blue";}
	if (X[i,5] == "Iris-virginica") {CC[i+1] = "green";}
}
V(g)$color = CC
V(g)$names = c(0,as.integer(X[,5]))
V(g)$clust = c(rep(1,vcount(g)))
V(g)$ids = seq(1:151)-1

IS = which(X[,5] == "Iris-setosa")

if (alg == "full"){
	for (i in 1:149){
		for (j in (i+1):150){
			if (M[i,j] > simthresh){
				g = add.edges(g, c(i,j), weight=M[i,j])
			}
		}
	}
} else if (alg == "knn"){
	for (i in 1:150){
		Oi = order(M[i,], decreasing=TRUE)
		for (j in 1:150){
			if (i != j && (j < knn || M[i,Oi[j]] > simthresh)){
				if (M[i,j] > simthreshKNN){
					g = add.edges(g, c(i,j), weight=M[i,j])
				}
			}
		}
	}
}

#delete all small components, and keep only the largest one
comps = clusters(g)
del_c = order(comps$csize, decreasing=TRUE)[3:comps$no]
g=delete.vertices(g,which(comps$membership %in% (del_c-1))-1)


if (file.exists('iris-base-layout.txt')){
	l = as.matrix(read.table('iris-base-layout.txt'))
	l[59,] = (l[10,]+l[60,])/2 #adjust the mismatched node
} else {
	l = layout.fruchterman.reingold(g, repulserad=200, niters=2000,
		area=500)
	l = layout.norm(l,0,5,0,5,0)
}

#make weighted adjacency matrix
#EE = get.edgelist(g)
#A = matrix(0,vcount(g),vcount(g))
#for (i in 1:ecount(g)){
	#A[EE[i,1]+1,EE[i,2]+1] = E(g)[i-1]$weight
	#A[EE[i,2]+1,EE[i,1]+1] = E(g)[i-1]$weight
#}

A = M[V(g)$ids,V(g)$ids]
A[A <= simthresh] = 0
#diag(A) = 0

K = 3
n = nrow(A)

#parition points into k clusters randomly
bs = ceiling(n/K)
SH <- sample(1:n)
W = matrix(0,n,0)
for(i in 1:K){
        Pi = c(rep(0, n))
        s = (i-1)*bs+1
		e = min(n,i*bs)
        Pi[SH[s:e]] = 1.0
		W = cbind(W,Pi)
}

#start Kmeans iterations
notconverged=TRUE
niters = as.integer(0)

oneK = matrix(1,K,1)
oneN = matrix(1,n,1)
while(notconverged){
	Wbn = t(t(W)/colSums(W))
	S = 2*(A %*% Wbn)
	T = oneN %*% diag(t(Wbn) %*% A %*% Wbn)
	D =  T-S
	AssignVec = apply(D,1,which.min)
	nW = matrix(0,n,0)
	for (i in 1:K){
		Pi = c(rep(0, n))
		Pi[which(AssignVec == i)] = 1.0
		nW = cbind(nW,Pi)
	}
	diff = sum(abs(W-nW))
	if (diff == 0 || niters > 50) {notconverged = FALSE}
	W = nW
	print(c("ITER", niters, diff))
	niters = niters + 1
}

colSums(W)
#compute SSE value
SSE = 0
for (i in 1:K){
	nz = which(W[,i] > 0)
	ni = length(nz)
	SSE = SSE + (sum(diag(A)[nz]) - sum(A[nz,nz])/ni)
}
SSE

#compute normalized cut value
Jw = 0
for (i in 1:K){
	Ci = which(AssignVec==i)
	voli = sum(A[Ci,])
	cuti = sum(A[Ci,-Ci])
	Jw = Jw + (cuti/voli)
}
print(c("normalized cut value", Jw))


#delete the cross-cluster edges
for (u in 1:vcount(g)){
	for (v in 1:vcount(g)){
		if (are.connected(g,(u-1),(v-1)) && AssignVec[u] !=
			AssignVec[v]){
			g = delete.edges(g,E(g,P=c(u-1,v-1)))
		}
	}
}

colary = c("red","blue","green")
V(g)$color = colary[AssignVec]

plot(g,layout=l, vertex.size=5, vertex.label=NA)

labs = c("Iris-setosa", "Iris-versicolor", "Iris-virginica")
fname = paste("iris-kernelkmeans-layout",".txt",sep="")
cat (file=fname)
for (v in V(g)){
    cat (file=fname, append=TRUE, "v", v, l[v+1,1],
            l[v+1,2], labs[V(g)$names[v+1]], 
            AssignVec[v+1], "\n")
}

EE = get.edgelist(g)
for (i in 1:ecount(g)){
    cat (file=fname,  append=TRUE, "e", EE[i,1], EE[i,2], "\n")
}

CK = matrix(0,4,K)
for (i in 1:K){
	Ci = which(AssignVec==i)
	CK[1,i] = length(which(labs[V(g)$names[Ci]]=="Iris-setosa"))
	CK[2,i] = length(which(labs[V(g)$names[Ci]]=="Iris-virginica"))
	CK[3,i] = length(which(labs[V(g)$names[Ci]]=="Iris-versicolor"))
	CK[4,i] = length(Ci)
}
print(c("normalized cut value, SSE", Jw, SSE))
print(CK)
