library(MASS)
library(mvtnorm)
library(gplots)

setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/probabilistic/figs')
X = read.table('iris.txt', sep=",")
D = as.matrix(X[1:2])

D1 = subset(D,X[,5] == 'Iris-setosa')
D2 = subset(D,X[,5] != 'Iris-setosa')
n1 = dim(D1)[1]
n2= dim(D2)[1]

m1 = colSums(D1)/n1
m2 = colSums(D2)/n2

S1 = (n1-1)/n1*cov(D1)
S2 = (n2-1)/n2*cov(D2)
S1
S2
dmvnorm(x=c(6.75,4.25),mean=m1,sigma=S1)
0.33*dmvnorm(x=c(6.75,4.25),mean=m1,sigma=S1)
dmvnorm(x=c(6.75,4.25),mean=m2,sigma=S2)
0.67*dmvnorm(x=c(6.75,4.25),mean=m2,sigma=S2)

dmvnorm(x=c(6.75,4.25),mean=m1,sigma=diag(diag(S1)))
0.33*dmvnorm(x=c(6.75,4.25),mean=m1,sigma=diag(diag(S1)))
dmvnorm(x=c(6.75,4.25),mean=m2,sigma=diag(diag(S2)))
0.67*dmvnorm(x=c(6.75,4.25),mean=m2,sigma=diag(diag(S2)))



#categorical
ejpmf = function(X,n){
	ba = c(4.2,5.2,6.1,7.0,7.9)
	bb = c(1.9,2.8,3.6,4.4)
	hA = hist(X[,1],breaks=ba,plot=F)
	hB = hist(X[,2],breaks=bb,plot=F)
	na = length(ba)
	nb = length(bb)

	ct = matrix(0,nrow=na-1, ncol=nb-1)
	for (i in 1:(na-1)){
		for (j in 1:(nb-1)){
			ct[i,j] = length(which(
					X[,1] > ba[i] & X[,1] <= ba[i+1] & 
					X[,2] > bb[j] & X[,2] <= bb[j+1]))
		}
	}
	show(ct)
	ct/n
}

ejpmf(D1,n1)
ejpmf(D2,n2)

