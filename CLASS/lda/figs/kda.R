library(methods)
library(MASS)
setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/lda/figs')

dataF = "iris"
dataF = "irisPC"

D1 = matrix(0,50, 2)
D2 = matrix(0,100, 2)
if (dataF == "iris"){
	X = read.table('iris.txt', sep=",")
	D1 = subset(X[1:2],X[,5] == 'Iris-virginica')
	D2 = subset(X[1:2],X[,5] != 'Iris-virginica')
} else if (dataF == "irisPC"){
	X = read.table('iris-PC.txt', sep=" ")
	D1 = subset(X[1:2],X[,3] != 1)
	D1
	D2 = subset(X[1:2],X[,3] == 1)
}

n1 = dim(D1)[1]
n2= dim(D2)[1]
D = as.matrix(rbind(D1,D2))

if (dataF == "iris"){
write.table(D, file="iris-virginica.dat", row.names=FALSE,
		col.names=FALSE)
}

range(D[,1])
range(D[,2]) 

n=n1+n2

#D
n

#kernel matrix
K = D%*%t(D)
K = K*K #homogeneous quadratic

diff = sum(K-t(K))
show(c("diff", diff))

#K

m1 = rowSums(K[,1:n1])/n1
m2 = rowSums(K[,(n1+1):n])/n2

M = (m1-m2)%*%t(m1-m2)

On1 = matrix(1/n1,n1,n1)
In1 = diag(n1)
On2 = matrix(1/n2,n2,n2)
In2 = diag(n2)

N1 = K[,1:n1] %*% (In1-On1) %*% t(K[,1:n1])
N2 = K[,(n1+1):n] %*% (In2-On2) %*% t(K[,(n1+1):n])

N = N1 + N2
dim(N)

NI = ginv(N)

E = NI %*% M
ee = eigen(E)
show("eigen")
ee$values[1:10]
l1 = ee$value[1]
l1
a= ee$vectors[,1]
sqrt(sum(a*a))

s = sqrt(t(a) %*% K %*% a)
s
a = a/s
a

#direct solution without eigenvectors
b = as.vector(NI %*% (m1-m2))
dim(b)
sb = sqrt(t(b) %*% K %*% b)
sb
b = b/sb
b

#full approach homogeneous
s2 = sqrt(2)
Z = matrix(0,n,3)
for (i in 1:n){
	x1 = D[i,1]
	x2 = D[i,2]
	Z[i,] = c(s2*x1*x2, x1*x1, x2*x2)
}

write.table(Z[1:n1,], file="iris-QHc1.dat", row.names=FALSE,
		col.names=FALSE)
write.table(Z[(n1+1):n,], file="iris-QHc2.dat", row.names=FALSE,
		col.names=FALSE)
show("range")
range(Z[,1])
range(Z[,2])
range(Z[,3])
mmZ = colSums(Z)/n
show("mean Z")
mmZ

u1 = colSums(Z*a)
u1
#u1/sqrt(sum(u1*u1))
uu1 = colSums(Z*b)
uu1

#project points on u1
dim(t(a))
a = as.matrix(as.real(a), 150,1)
P = K %*% a
range(P)
show("direction vector range")
u1*(range(P)[1]*1.1)
u1*(range(P)[2]*1.1)

OO = matrix(0,n,1)
PROJ = cbind(P,OO)
#PROJ
write.table(PROJ, file="projw.dat", row.names=FALSE,
		col.names=FALSE)
pm1 = mean(PROJ[1:n1])
pm1*u1
ps1 = var(PROJ[1:n1])*(n1-1)
pm2 = mean(PROJ[(1+n1):n])
pm2*u1
ps2 = var(PROJ[(1+n1):n])*(n2-1)
show("m1, m2, s1, s2")
pm1
pm2
ps1
ps2
Jw = (pm1-pm2)^2/ (ps1+ps2)
show(c("Jw", Jw))

meanZ = t(matrix(rep(as.real(mmZ),n), 3, n))
UU = t(matrix(rep(as.real(u1),n), 3, n))
show("map")
dim(P)
MAP = UU*c(P)
#MAP
write.table(MAP[1:n1,], file="kdamappedPoints-c1.dat", 
		row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(MAP[(n1+1):n,], file="kdamappedPoints-c2.dat", 
		row.names=FALSE, col.names=FALSE, quote=FALSE)

LINES = c(rep('',n))
for (i in 1:n){
	zp = paste(round(Z[i,],digits=3),collapse=",")
	mp = paste(round(MAP[i,],digits=3),collapse=",")
	LINES[i] = paste("\\pstThreeDLine[linecolor=gray]",
			"(", zp, ")", "(", mp, ")",sep="",collapse="")
}
write.table(LINES, file="kdaproject.tex", row.names=FALSE,
		col.names=FALSE, quote=FALSE)

#call inbuilt LDA for verification
#Z
L = c(rep('c1',n))
L[(n1+1):n] = 'c2'
#L
ll = lda(Z,L)
ll
u1 = ll$scaling
u1/sqrt(sum(u1*u1))

