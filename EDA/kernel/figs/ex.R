#!/usr/bin/env Rscript
library(MASS)

M = t(matrix(c(
  5.9, 3, 
  6.9, 3.1,
  6.6, 2.9,
  4.6, 3.2,
  6, 2.2), 2, 5))

print(M)
mu = colSums(M)/5
print(mu)
print(sum(mu*mu))

LK = M%*%t(M)
QK = LK*LK

print("kernels: linear and quadratic")
print(LK)
print(sum(LK))
print(mean(LK))
print(QK)

#tot var
tv = 1/5*sum(diag(LK)) - mean(LK)
print("total var")
print(tv)

#centered linear
oo = matrix(1/5,5,5)
I = diag(5)
zz = I-oo
print(zz)
zLK = zz %*% LK %*% zz
print ("centered linear")
print(zLK)

Z = scale(M, center=TRUE, scale=FALSE)
print(Z)
zLK = Z%*%t(Z)
print(zLK)

#normalized
Dq = matrix(diag(QK),5,1)
print(Dq)
ss = sqrt(Dq %*% t(Dq))
print(sqrt(ss))
nQK = QK/ss
print("normalized Q")
print(nQK)

D = matrix(diag(LK),5,1)
print(D)
ss = sqrt(D %*% t(D))
print(sqrt(ss))
nLK = LK/ss
print("normalized")
print(nLK)

#centered and normalized
print("normalized and centered")
Ws = diag(1/sqrt(diag(zLK)))
ncLK = Ws %*% zLK %*% Ws
print(ncLK)

#scale/unitize the LK matrix
ss = sqrt(rowSums(M*M))
sM = scale(t(M), center=FALSE, scale=ss)
sM = t(sM)
print(sM)
sLK = sM%*%t(sM)
print(sLK)

#norm, distance
d12 = QK[1,1] + QK[2,2] - 2*QK[1,2]
print(d12)
print(sqrt(d12))

#gaussian kernel
GK = matrix(0,5,5)
for (i in 1:5){
  for (j in 1:5){
	dd = M[i,]-M[j,]
	GK[i,j] = exp(-1*sum(dd*dd)/2)
  }
}
print("Gaussian")
print(GK)
print(sum(GK))
mGK = mean(GK)
print(mGK)
print(sqrt(sum(GK))/5)

d1mGK = GK[1,1] - 2/5*sum(GK[1,]) + mGK
print(d1mGK)
print(sum(GK[1,]))

#eigen
print("eigendecomp LK")
eLK  = eigen(LK)
print(eLK)
U = eLK$vectors
L = sqrt(eLK$values[1:2])
print(L)
Px1 = (L[1]*U[,1])
Px2 = (L[2]*U[,2])
print(Px1)
print(Px2)
print(c(Px1[1]*Px1[2], Px2[1]*Px2[2]))
print(Px1[1]*Px1[2] + Px2[1]*Px2[2])

print("inverse matrices for LK")
ILK = ginv(LK)
print(ILK)
dd = diag(1/L^2)
print(dd)
iiLK = U[,1:2] %*% dd %*% t(U[,1:2])
print(iiLK)

sILK = U[,1:2] %*% diag(1/L) %*% t(U[,1:2])
print("empirical map")
Emap = (sILK) %*% LK
print(Emap)
print(Emap%*%Emap)

print("dot prod of rows of LK")
LKd = LK %*% LK
print(LKd)
ee = eigen(LKd)
print(ee)
