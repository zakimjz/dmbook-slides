setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/svm/figs')
X <- read.table("iris.txt", sep=",")
D <- X[1:2]
l <- c(rep(1,150))
l[which(X[,5] == "Iris-setosa")] = -1
D = cbind(D,l)
D1 = D[(which(X[,5] == "Iris-setosa")),]
D2 = D[(which(X[,5] != "Iris-setosa")),]
Dx = rbind(D1,D2)
write.table(Dx, "iris-slwc.txt", sep=" ",row.names=F,col.names=F)

pr <- prcomp(X[1:4])
l <- c(rep(1,150))
l[which(X[,5] == "Iris-versicolor")] = -1
Z <- cbind(pr$x[,1:2],l)
D1 = Z[(which(X[,5] == "Iris-versicolor")),]
D2 = Z[(which(X[,5] != "Iris-versicolor")),]
Dx = rbind(D1,D2)
write.table(Dx, "iris-PC.txt", sep=" ",row.names=F,col.names=F)
