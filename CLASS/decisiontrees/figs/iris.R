library(methods)
setwd('/Users/zaki/research/DataMiningBook/dm08/CLASS/decisiontrees/figs')


entropy <- function(nP,nN){
	pP = nP/(nP+nN)
	pN = nN/(nP+nN)
	tP = 0
	tN = 0
	if (pP > 0){ tP = pP*log2(pP) }
	if (pN > 0){ tN = pN*log2(pN) }
	return (-1*(tP+tN))
}

splitEntropy <- function(nPL,nNL,nPR,nNR){
	nL = nPL+nNL
	nR = nPR+nNR
	n = nL+nR
	eL = eR = 0
	if (nL > 0){ eL  = (nL/n)*entropy(nPL,nNL) }
	if (nR > 0){ eR  = (nR/n)*entropy(nPR,nNR) }

	return (eL+eR)
}


bestSplit <- function(AV,level){
	n = dim(AV)[1]
	U = unique(AV[,2])
	nU = length(U)
	U = (U[1:(nU-1)]+U[2:nU])/2
	minE = 1
	minV = -1
	nP = length(which(AV[,3] == 1))
	nN = length(which(AV[,3] == -1))
	E = entropy(nP,nN)
	if (level == 0){
		G = vector()
		P = vector()
		N = vector()
	}
	for (v in U){
		nPL = length(which(AV[,2] <= v & AV[,3] == 1))
		nNL = length(which(AV[,2] <= v & AV[,3] == -1))
		nPR = length(which(AV[,2] > v & AV[,3] == 1))
		nNR = length(which(AV[,2] > v & AV[,3] == -1))
		splitE = splitEntropy(nPL,nNL,nPR,nNR)
		gain = E - splitE
		#print(c(v,nPL,nNL,nPR,nNR,splitE))
		if (level == 0){
			G = append(G,gain)	
			P = append(P,nPL)
			N = append(N,nNL)
		}
		if (splitE < minE){
			minE = splitE
			minV = v
		}
	}
	if (level == 0){#print the gain values and cumulative dist
		UG = cbind(U, G, P, N)
		print(UG)
		print(c(nP, nN))
	}
	return (list(V=minV,E=minE))
}

DT <- function(AV1,AV2,level){
	#print(AV1)
	#print(AV2)
	n = length(AV1)/3
	#print(n)
	if (n == 1){
		AV1 = t(as.matrix(AV1))
		AV2 = t(as.matrix(AV2))
	}
	nP = length(which(AV1[,3] == 1))
	nN = length(which(AV1[,3] == -1))
	purity = max(nP/n, nN/n)
	indent = c(rep("\t",level))
	cat(c(indent, level, n, nP, nN, purity))
	if (n > minsize & purity < minpurity){
		bs1 = bestSplit(AV1,level)
		bs2 = bestSplit(AV2,level)
		#print(bs1)
		#print(bs2)
		if (bs1$E < bs2$E){ 
			bs = bs1 
			Lidx = AV1[which(AV1[,2] <= bs$V),1]
			Ridx = AV1[which(AV1[,2] > bs$V),1]
			attr ="AV1"
		} else {
			bs = bs2 
			Lidx = AV2[which(AV2[,2] <= bs$V),1]
			Ridx = AV2[which(AV2[,2] > bs$V),1]
			attr = "AV2"
		}
		show(c(attr, bs$E, bs$V))
		
		level = level+1
		#split the DT AV lists into two
		AV1L = AV1[AV1[,1] %in% Lidx,]		
		AV1R = AV1[AV1[,1] %in% Ridx,]		
		AV2L = AV2[AV2[,1] %in% Lidx,]		
		AV2R = AV2[AV2[,1] %in% Ridx,]		
		DT(AV1L, AV2L, level)
		DT(AV1R, AV2R, level)
	} else {
		cat("\n")
	}
}


X <- as.matrix(read.table('iris-slwc.txt'))
n <- dim(X)[1]
d <- dim(X)[2]
Y <- X[,3]

ID <- seq(1:n)
X1 <- cbind(ID,X[,1],Y)
X2 <- cbind(ID,X[,2],Y)
X1 <- X1[order(X1[,2]),]
X2 <- X2[order(X2[,2]),]

minsize = 1
minpurity = 0.95
level = 0
DT(X1,X2,level)

#compute the gain for categorical splits

pmf <- function(X1, V){
	nPL = length(which(X1[,2] %in% V & X1[,3]==1))
	nNL = length(which(X1[,2] %in% V & X1[,3]==-1))
	nPR = length(which(!(X1[,2] %in% V) & X1[,3]==1))
	nNR = length(which(!(X1[,2] %in% V) & X1[,3]==-1))
	return(c(nPL, nNL, nPR, nNR))
}
X1[which(X1[,2] >= 4.3 & X1[,2] <= 5.2),2] = "a1"
X1[which(X1[,2] > 5.2 & X1[,2] <= 6.1),2] = "a2"
X1[which(X1[,2] > 6.1 & X1[,2] <= 7.0),2] = "a3"
X1[which(X1[,2] > 7.0 & X1[,2] <= 7.9),2] = "a4"

f = pmf(X1, c('a1','a2','a3','a4'))
eD = splitEntropy(f[1],f[2],f[3],f[4])
print(c(f, eD))

f = pmf(X1, c('a1'))
eS = splitEntropy(f[1],f[2],f[3],f[4])
print(c(f,eS, eD-eS))

f = pmf(X1, c('a2'))
eS = splitEntropy(f[1],f[2],f[3],f[4])
print(c(f,eS, eD-eS))

f = pmf(X1, c('a3'))
eS = splitEntropy(f[1],f[2],f[3],f[4])
print(c(f,eS, eD-eS))

f = pmf(X1, c('a4'))
eS = splitEntropy(f[1],f[2],f[3],f[4])
print(c(f,eS, eD-eS))

f = pmf(X1, c('a1', 'a2'))
eS = splitEntropy(f[1],f[2],f[3],f[4])
print(c(f,eS, eD-eS))

f = pmf(X1, c('a1', 'a3'))
eS = splitEntropy(f[1],f[2],f[3],f[4])
print(c(f,eS, eD-eS))

f = pmf(X1, c('a1', 'a4'))
eS = splitEntropy(f[1],f[2],f[3],f[4])
print(c(f,eS, eD-eS))

f = pmf(X1, c('a2', 'a3'))
eS = splitEntropy(f[1],f[2],f[3],f[4])
print(c(f,eS, eD-eS))

f = pmf(X1, c('a2', 'a4'))
eS = splitEntropy(f[1],f[2],f[3],f[4])
print(c(f,eS, eD-eS))

f = pmf(X1, c('a3', 'a4'))
eS = splitEntropy(f[1],f[2],f[3],f[4])
print(c(f,eS, eD-eS))

