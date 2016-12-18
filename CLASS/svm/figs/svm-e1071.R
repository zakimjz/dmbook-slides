#!/usr/bin/env Rscript
library(e1071)

h <- function(x){
	return (t(w) %*% x - b)
}


#homogeneous quadratic phi mapping
phi2h <- function(x){
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


args <- commandArgs(TRUE)
#print(args)
datafile = args[1]
C <- as.integer(args[2]) #regularization constant

kernel = 'linear'
if (length(args) > 2){
        kernel = args[3]
}


D <- as.matrix(read.table(datafile))
n <- nrow(D)
d <- ncol(D)

X <- D[,1:(d-1)]
Y <- D[,d]
ones <- c(rep(1,nrow(D)))
X <- as.matrix(cbind(X,ones))

if (kernel == 'linear'){
	model <- svm(X, Y, type="C-classification",
	kernel="linear", cost=C, scale=F);
} else if (kernel == 'quadratic'){
	model <- svm(X, Y, type="C-classification",
	kernel="polynomial", degree=2, cost=C, gamma=1, coef0 = 0, scale=F)
}

#summary(model)
#print(model$SV)
#print(model$rho)

print(model$index)
print(model$coefs)
AY = c(rep(0,n)) 
AY[model$index] = model$coefs
if (kernel == 'linear'){
	w <- t(X) %*% AY #compute the weight vector
} else if (kernel == 'quadratic'){
    phiX = apply(X, 1, phi2h)
    w = phiX %*% AY
}

#w <- t(model$SV) %*% model$coefs
b <- model$rho


cat("Weight Vector\n")
print (w)
cat("\nBias\n")
print (b)

#preds <- apply(X, 1, function(x) return (ifelse(h(x)>0, 1, -1)))
#pp <- preds*Y
#errors <- which(pp < 0)
#cat ("\nError Rate:", length(errors)/n, "\n")

pred <- as.vector(as.integer(predict(model,X)))
pred[which(pred==1)]=-1
pred[which(pred==2)]=1
pp <- pred*Y
errors <- which(pp < 0)
cat ("\nError Rate:", length(errors)/n, "\n")

