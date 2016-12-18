#!/usr/bin/env Rscript

#homogeneous quadratic phi mapping
phi2h <- function(x,d){
  r = c(rep(0,d*(d+1)/2));
  k=1;
  for (i in 1:d){
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

############MAIN##########
setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/eval/figs')
D = read.table('iris-PC.txt')
n = dim(D)[1]
d = dim(D)[2]
X = D[,1:(d-1)]
phiX = apply(X, 1, phi2h, 2)
phiX = t(phiX)
#print(phiX)
print(apply(phiX, 2, range)) 
#write.table(phiX, "iris-PC-3D.txt", row.names=FALSE, col.names=FALSE)
for (i in 1:n){
  str = sprintf("\\psPoint(%0.2f, %0.2f, %0.2f){p%d}", 
				phiX[i,1], phiX[i,2], phiX[i,3], i)
  str2 = sprintf("\\psdot(p%d)", i)
  print(str)
  print(str2)
}
