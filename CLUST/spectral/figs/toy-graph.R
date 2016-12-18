setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/spectral/figs')
library(igraph)
library(Cairo)

A = matrix(c(0,1,0,1,0,1,0, 
			 1,0,1,1,0,0,0, 
			 0,1,0,1,0,0,1, 
			 1,1,1,0,1,0,0, 
			 0,0,0,1,0,1,1, 
			 1,0,0,0,1,0,1, 
			 0,0,1,0,1,1,0), 7,7)
n=7
I = diag(n)


#RUN MCL
Am = A+I
ID = diag(1/rowSums(Am))
M = ID%*%Am
print(format(Am,digits=2))
cs = rowSums(M)
Z = M/cs
oZ = Z;
eps = 0.001
p = 2.5 #powers of the elements
iter = as.integer(1)
diff = 1
while(diff > eps){
	print(c("iter", iter))
	print(format(Z,digits=2))
	Z = Z%*%Z
	print(c("Z*Z"))
	print(format(Z,digits=2))
	Z = Z^p
	cs = rowSums(Z)
	Z = Z/cs
	print(c("Z^r"))
	print(format(Z,digits=2))
	diff = sum((Z-oZ)^2)
	print(c(iter, diff))
	oZ = Z
	iter = iter+1
}

Z[which(Z < 10^{-5})] = 0 #make very small values 0
Z[which(Z > 0)]=1

ng = graph.adjacency(Z,mode="directed",weighted=TRUE)
clust = clusters(ng)
clust

#NORM AND RATIO CUT
print(format(A,digits=2))
laplace = "normalized"
#laplace = "ratio"
D = diag(rowSums(A))

#asymmetric laplacian
ID = diag(diag(D^{-1}))
L = ID %*% (D-A)
print(format(L,digits=2))

#modularity
d = as.matrix(rowSums(A),n,1)
sd = sum(d)
E = d%*%t(d)
E = E/(sd^2)
L = A/sd-E


if (laplace == "ratio"){
	L = D-A
	#l[69,2]=l[91,2] #adjust layout for one of the cluster nodes
	#l[126,]=c(1.7,2.2)
} else {
	ID = diag(diag(D^{-1/2}))
	L = ID %*% (D-A) %*% ID
}


#get the smallest eigenvectors
ss=eigen(L)
n = 7
k = 2
U = ss$vectors[,n:(n-k+1)]
print(format(U,digits=3))

#U = ss$u[,1:k]
for (i in 1:n){
	U[i,] = U[i,]/sqrt(sum(U[i,]*U[i,]))
}
print(format(U,digits=3))

kk = kmeans(U,k)



