setwd('/Users/zaki/research/DataMiningBook/dm08/CLUST/representative/figs')

library(mvtnorm)

#two normals
m1=c(3,5)
#s1=matrix(c(1,0.4,0.4,0.2),2,2)
s1=matrix(c(1,0,0,0.1),2,2)

m2=c(10,3)
s2=matrix(c(1,-0.3,-0.3,0.2),2,2)

C1 = rmvnorm(100,m1,s1)
C2 = rmvnorm(100,m2,s2)

# one parabolic: tip at (h,k), width by a
a=0.15
h=3
k=2.5
C3 = matrix(0,100,2)
for (i in 1:100){
	x = runif(1,0.25,7)
	y = a*(x-h)^2+k
	x = x + runif(1,0,0.5) #add noise
	y = y + runif(1,0,0.5) #add noise
	C3[i,] = c(x,y) #scale the y coord
}

D = rbind(C1,C2,C3)
plot(D)

Y = c(rep(1,100), rep(2,100), rep(3,100))
Z = cbind(D,Y)
write.table(Z,'nldata.txt', row.names=FALSE, col.names=F, sep=',')
