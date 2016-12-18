library(gplots)
setwd('/Users/zaki/research/DataMiningBook/dm08/EDA/categorical/figs')
X <- read.table('iris.txt', sep=",")

variances <- function(CC,n){
	nn = length(CC)
	ones = c(rep(1,nn))
	m = CC/n
	show(m)
	mm = ones-m
	show(mm)
	M = diag(m)
	U = m %*% t(m)
	show(M-U)
	M-U
}

h2d <- function(A,B){
	m = length(A)
	cA = A
	cB = B
	for (i in 1:m){
		if (A[i] >= 4.3 && A[i] <= 5.2){cA[i] = "a11"}
		else if (A[i] > 5.2 && A[i] <= 6.1){cA[i] = "a12"}
		else if (A[i] > 6.1 && A[i] <= 7.0){cA[i] = "a13"}
		else if (A[i] > 7.0){cA[i] = "a14"}
		if (B[i] >= 2.0 && B[i] <= 2.8){cB[i] = "a21"}
		else if (B[i] > 2.8 && B[i] <= 3.6){cB[i] = "a22"}
		else if (B[i] > 3.6){cB[i] = "a23"}
	}
	ct = table(cA,cB)
	ct
}

covar <- function(A,B, CC1, CC2, n){
	nn1 = length(CC1)
	m1 = CC1/n
	nn2 = length(CC2)
	m2 = CC2/n
	M12  = m1 %*% t(m2)

	ct = matrix(0,nn1,nn2)
	#hAB = hist2d(A,B,nbins=c(nn1,nn2),show=F)
	hAB = h2d(A,B)
	#ct = as.matrix(hAB$counts)
	ct[1:nrow(hAB),1:ncol(hAB)] = hAB
	show(ct)
	show(ct/n)

	show(ct/n-M12)
	ct
}

analyze <- function(X,n){
	A <- X[,1]
	B <- X[,2]
	C <- X[,6]
	hA = hist(A,breaks=(c(4.3,5.2,6.1,7.0,7.9)),plot=F)
	hB = hist(B,breaks=(c(2,2.8,3.6,4.4)),plot=F)
	hC = hist(C,breaks=(c(0,1,2,3)),plot=F)

	variances(hA$counts, n)
	variances(hB$counts, n)
	variances(hC$counts, n)

	covar(A,B,hA$counts,hB$counts,n)
	covar(A,C,hA$counts,hC$counts,n)
	covar(B,C,hB$counts,hC$counts,n)

	#contingency for SL and SW
	ct = covar(A,B,hA$counts,hB$counts,n)

	rs = rowSums(ct)
	cs = colSums(ct)
	show(rs)
	show(cs)
	show(rs%*%t(cs))
	E = rs%*%t(cs)/n
	show(E)
	XS = (ct-E)^2/E
	sum(XS)
}

countmatrix <- function (A,B,n){
	hA = hist(A,breaks=(c(4.3,5.2,6.1,7.0,7.9)),plot=F)
	hB = hist(B,breaks=(c(2,2.8,3.6,4.4)),plot=F)
	ct = covar(A,B,hA$counts,hB$counts,n)
	ct
}

ve = as.numeric(subset(rownames(X),X$V5=="Iris-versicolor"))
se = as.numeric(subset(rownames(X),X$V5=="Iris-setosa"))
vi = as.numeric(subset(rownames(X),X$V5=="Iris-virginica"))
Z = as.vector(c(rep(0,150)))
Z[ve] = 1
Z[se] = 2
Z[vi] = 3
X <- cbind(X,Z)

#3D analysis
analyze(X,150)

#for class
m2 = c(rep(1/3,3))
M2 = diag(m2)
U2 = m2 %*% t(m2)
show(M2-U2)

#3-way contingency
SE = X[se,]
VE = X[ve,]
VI = X[vi,]
c1 = countmatrix(VE[,1],VE[,2],150)
c2 = countmatrix(SE[,1],SE[,2],150)
c3 = countmatrix(VI[,1],VI[,2],150)

A <- X[,1]
B <- X[,2]
C <- X[,6]
hA = hist(A,breaks=(c(4.3,5.2,6.1,7.0,7.9)),plot=F)
hB = hist(B,breaks=(c(2,2.8,3.6,4.4)),plot=F)
hC = hist(C,breaks=(c(0,1,2,3)),plot=F)

Eab = (hA$counts %*% t(hB$counts))*hC$counts[1]/(150^2)

NN = cbind(c1, c2, c3)
EE = cbind(Eab, Eab, Eab)
XS = (NN-EE)^2/EE
sum(XS)
