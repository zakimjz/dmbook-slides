setwd('/Users/zaki/research/DataMiningBook/dm08/EDA/graph/figs')
library(igraph)
library(Cairo)
X = read.table('iris.txt', sep=",")
D = as.matrix(X[1:4])
M = matrix(0,150,150);
#for (i in 1:149){
	#for (j in (i+1):150){
		#dd = D[i,]-D[j,];
		#sdd = sum(dd*dd);
		#dij  = exp(-1*sdd);
		#M[i,j] = dij;
	#}
#}

for (i in 1:150){
	Zi = t(t(D) - D[i,]);
	rsi = rowSums(Zi*Zi);
	di = exp(-1*rsi);
	M[i,] = di;
}

Z = M[M>0]
m = mean(Z)
s= sd(Z)

simthresh = m + 2*s

g = graph.empty(n=150, directed=F)

for (i in 1:149){
	for (j in (i+1):150){
		if (M[i,j] > simthresh){
			cat ("e", i, j, M[i,j], "\n")
			g = add.edges(g, c(i,j))
		}
	}
}

l = layout.fruchterman.reingold(g, repulserad=5000, niters=2000,
		area=3000)
l = layout.norm(l,0,5,0,5,0)

labs = c("Iris-setosa", "Iris-versicolor", "Iris-virginica")
for (i in 1:150){
	cat ("v", i, l[i,1], l[i,2], labs[X[i,5]], "\n")
}

CC = c(rep(0,150))
for (i in 1:150){
	if (X[i,5] == "Iris-setosa") {CC[i] = "red";}
	if (X[i,5] == "Iris-versicolor") {CC[i] = "blue";}
	if (X[i,5] == "Iris-virginica") {CC[i] = "green";}
}

V(g)$color = CC

d = get.diameter(g)
E(g)$color <- "grey"
E(g)$width <- 1
E(g, path=d)$color <- "red"
E(g, path=d)$width <- 2
#plot(g,layout=l, vertex.size=5, vertex.label=NA)

cat("d", d, "\n")

#degree.distribution(g)
dd =degree(g)
degdist = c(rep(0,max(dd)+1))
for (i in 1:length(dd)){
  degdist[dd[i]+1] = degdist[dd[i]+1]+1
}
sum(degdist)
print(degdist)
pp = degdist/sum(degdist)
pp
for (i in 1:(max(dd)+1)){
  print(c(as.integer(i-1), pp[i]))
}
ph = path.length.hist(g)
ph$res 
sum(ph$res) 
cumsum(ph$res)/sum(ph$res)
mm = seq(1,11)
sum(ph$res*mm)
sum(ph$res*mm)/sum(ph$res)
