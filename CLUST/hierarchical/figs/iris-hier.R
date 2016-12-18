setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/hierarchical/figs')
library(igraph)
library(Cairo)

####MAIN####
k = 3 #num clusters
method="single"
#"ward", "single", "complete", "average", "mcquitty", 
# "median" or "centroid".

args <- commandArgs(trailingOnly = TRUE)
k = as.integer(args[1])
method = args[2]

print(c(k, method))
X = read.table('iris-PC.txt', sep=",")
D = as.matrix(X[1:2])
n = nrow(D);
classattr = 3
dd <- dist(D, method = "euclidean")
fit = hclust(dd, method=method)
plot(fit, labels=X[,classattr])
kk <- cutree(fit, k=k)
print(kk)

#print(fit)
#print(fit$merge)
print(table(X[,classattr], kk))

YY = cbind(D,kk)
#print(YY)
labs  = c("Iris-setosa", "Iris-versicolor", "Iris-virginica")
ofname = 'results.txt'
idx =0
write.table('%', ofname, row.names=FALSE, col.names=FALSE, quote=FALSE)
for (i in 1:3){
  for (j in 1:3){
	cC  = which(kk == i)
	cT = which(X[,classattr] == labs[j])
	iCT = intersect(cC,cT)
	if (length(iCT) > 0){
	  print (c(i, labs[j], length(iCT), idx))
	  ll  = '%'
	  ll = paste(ll, c(i, labs[j], length(iCT), idx))
	  #print(ll)
	  write.table(ll, ofname, 
				  row.names=FALSE, col.names=FALSE,
				  append=TRUE,quote=FALSE)
	  write.table(YY[iCT,], ofname, row.names=FALSE,
				col.names=FALSE, sep = " ", append=TRUE,quote=FALSE)
	  idx = idx + length(iCT)
	}
  }
}


