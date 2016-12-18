setwd('/Users/zaki/research/DataMiningBook/dm08/EDA/data/figs')
X = read.table('iris.txt', sep=",")
D = as.matrix(X[1:2])

m = colMeans(D)
m

o = c(rep(1,150))
Z = D-o%*%t(m)
Z

write.table(Z,"iris-slw-c.dat", row.names=F, col.names=F)

dd = D %*% t(D)

s = mean(diag(dd)) - t(m) %*% m
s

library(MASS)
ll = lda(Z)
ll
u = ll$scaling[,1]
mg = sqrt(t(u) %*% u)
u/mg

rownames(Z) = X[,5]
write.table(Z,"iris-slw-clbl.txt", col.names=F, quote=F)
