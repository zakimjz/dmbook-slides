#!/usr/bin/env Rscript
A = matrix(
	c(0, 0, 1, 1, 0,
	0, 0, 1, 0, 1, 
	1, 1, 0, 1, 0,
	1, 0, 1, 0, 1,
	0, 1, 0, 1, 0), 5,5)
A = t(A)

D = diag(rowSums(A))

nL = A-D
print(nL)

E = eigen(nL)
print(E)

L = E$values
U = E$vectors

K = matrix(0,5,5)
beta = 0.2
for (i in 1:5){
  #print(exp(beta*L[i]))
  #print(t(t(U[,i])) %*% t(U[,i]))
  K = K + exp(beta*L[i]) * (t(t(U[,i])) %*% t(U[,i]))
}
print(K)

EE = eigen(K)
print(EE)

#von Neumann Kernel
I = diag(5)
beta = 0.2
M = diag(L)
vS = I - beta*M
print(vS)
IvS = solve(vS)
print(IvS)
vK = solve(I-beta*nL)
print(vK)

#verify
vKK = U%*%IvS%*%t(U)
print(vKK)

#what about Laplacian instead of neg Laplacian
#ivKL = solve(I+beta*nL)
#print(vKL)
#print(eigen(vKL))


#M = diag(c(0,-2,-4,-5,-20))
#AA = U %*% M %*% t(U)
#beta = -0.049
#vS = I - beta*M
#print(vS)
#print(solve(vS))
#vK = solve(I-beta*nL)
#print(vK)

