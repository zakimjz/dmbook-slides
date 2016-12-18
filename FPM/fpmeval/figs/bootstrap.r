#!/usr/bin/env Rscript
library(arules)
setwd('/Users/zaki/research/DataMiningBook/dm08/FPM/fpmeval/figs')

confinterval = function(A){
  lb = quantile(A,(1-alpha)/2,type=1)
  ub = quantile(A,(1+alpha)/2,type=1)
  return (c(lb,ub))
}

pval = function(v,A){
  return (ecdf(A)(v))
}


#MAIN#
set.seed(2)
X = read.table('iris-discrete.txt', sep=',')
D = as(X, "transactions")
N = 150
K = 100
minsup=0.1
alpha = 0.9

#frequent itemsets in D
F = eclat(D, parameter = list(supp = minsup))
S = support(F, D)

#bootstrap samples
SS = c()
for (i in 1:K){
  Ds = sample(D, N, replace=TRUE)
  Ss = support(F, Ds)
  SS = cbind(SS, Ss)
}

SSS = cbind(S,SS)
#M = sapply(SS, mean)
#M

minpv = 1.0
mini = 0
for (i in 1:length(F)){
  minv = min(SS[i,])
  maxv = max(SS[i,])
  ci = confinterval(SS[i,])
  pv = pval(S[i], SS[i,])
  spv = pval(minsup, SS[i,])
  iset = as(items(F[i]),"list")
  if (pv < minpv){
	minpv = pv
	mini = i
  }
  #print(iset)
  #cat(i,S[i], ci, pv, spv, minv, maxv, "\n")
}
print (as(items(F[mini]),"list"))
cat("min", mini, S[mini], confinterval(SS[mini,]), pval(S[mini], SS[mini,]),
	pval(minsup,SS[mini,]), mean(SS[mini,]), "\n")

#EPMF
V = sort(unique(SS[mini,]))
VC = table(SS[mini,])/K
VPMF = cbind(V,as.matrix(VC))
VPMF
write.table(VPMF, file="epmf-bootstrap.txt", row.names=FALSE, 
			col.names=FALSE)

#ECDF
E = ecdf(SS[mini,])
PV = E(V)
VPV = cbind(V,PV)
print(VPV)
write.table(VPV, file="ecdf-bootstrap.txt", row.names=FALSE, 
			col.names=FALSE)
