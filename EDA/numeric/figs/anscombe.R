X1 = c(12,14,18,23,27,28,34,37,39,40)
X2 = c(3,5,10,20,35,40,43,60,25,27)*100
D = cbind(X1,X2)

sqrt(sum((D[1,]-D[2,])^2))

#range norm
S1 =sapply(X1, FUN=function(x) (x-min(X1))/(max(X1)-min(X1)))
S2 =sapply(X2, FUN=function(x) (x-min(X2))/(max(X2)-min(X2)))

S=cbind(S1,S2)
cov(S)*9/10
cor(S)

format(S1,digits=2)
format(S2,digits=2)

#z-norm
s1 = sqrt(var(X1)*9/10)
s2 = sqrt(var(X2)*9/10)
m1 = mean(X1)
m2 = mean(X2)

S1 =sapply(X1, FUN=function(x) (x-m1)/s1)
S2 =sapply(X2, FUN=function(x) (x-m2)/s2)

format(S1,digits=2)
format(S2,digits=2)

S=cbind(S1,S2)
cov(S)*9/10
cor(S)

