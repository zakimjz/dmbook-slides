library(MASS)
setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/lda/figs')
X = read.table('iris.txt', sep=",")
D = as.matrix(X[1:2])

mean(D[,1])
mean(D[,2])
D1 = subset(D,X[,5] == 'Iris-setosa')
D2 = subset(D,X[,5] != 'Iris-setosa')
n1 = dim(D1)[1]
n2= dim(D2)[1]
n = n1+n2

m1 = colSums(D1)/n1
m2 = colSums(D2)/n2

n1 
n2
m1
m2

m1-m2

B = (m1-m2)%*%t(m1-m2)
B

S1 = (n1-1)*cov(D1)
S2 = (n2-1)*cov(D2)
S = S1+S2
S
SI = ginv(S)
SI
M = SI%*%B
M
eigen(M)

w = SI%*%(m1-m2)
w
sqrt(sum(w*w))
w = w/sqrt(sum(w*w))
w

rownames(D) = X[,5]
write.table(D,"iris-slw-clbl.txt", col.names=F, quote=F)

L = c(rep('c1',n))
L[which(X[,5] != 'Iris-setosa')] = 'c2'
#L
ll = lda(D,L)
ll
u1 = ll$scaling
u1/sqrt(sum(u1*u1))

