assignpt <- function(X,m1,m2){
	c1 <- c()
	c2 <- c()
	a <- 1
	b <- 1
	for (i in 1:length(X)){
		if (abs(m1-X[i]) <= abs(m2-X[i])){
			c1[[a]] = X[i]
			a = a+1
		}
		else{
			c2[[b]] = X[i]
			b = b+1
		}
	}
	print(c1)
	print(c2)
	newm1 <- mean(c1)
	newm2 <- mean(c2)
	print(c(newm1,newm2))
	return (list(nm1=newm1, nm2=newm2))
}

X <- c(2,3,4,10,11,12,20,25,30)
m1 <- 2
m2 <- 4
nm = assignpt(X,m1,m2)
nm = assignpt(X,nm$nm1, nm$nm2)
nm = assignpt(X,nm$nm1, nm$nm2)
nm = assignpt(X,nm$nm1, nm$nm2)
nm = assignpt(X,nm$nm1, nm$nm2)

