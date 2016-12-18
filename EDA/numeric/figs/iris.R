setwd('/Users/zaki/research/DataMiningBook/dm08/EDA/numeric/figs')
D <- read.table('iris.txt', sep=',')
X <- D[,-5]
m = colSums(X)/150
m

Z = t(t(X)-m)

S = t(Z) %*% Z
S/150

cov(Z)
cov(X)

S12 = S[1:2,1:2]/150

eigen(S12)

W = eigen(S12)$vectors
L = diag(eigen(S12)$values)
u1 = W[,1]
u2 = W[,2]
e1 = c(1,0)
a= sum(u1*e1)
acos(a)*180/pi

#for plotting the u1 axis
s = (2-m[1])/u1[1]
u1[2]*s+m[2]

s = (9-m[1])/u1[1]
u1[2]*s+m[2]

#for u2
s = (1-m[2])/u2[2]
u2[1]*s+m[1]

s = (5-m[2])/u2[2]
u2[1]*s+m[1]

