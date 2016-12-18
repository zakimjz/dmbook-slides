library(mvtnorm)
d=100
m = rep(0,d)
n=100000

X = rmvnorm(n,m)

r = sqrt(d)
l = r-1/sqrt(2)
u = r+1/sqrt(2)

cnt = 0
for (i in 1:n){
	mag = sqrt(sum(X[i,]^2))
	if (mag >= l && mag <= u){
		cnt = cnt+1
	}
}

print(cnt)
