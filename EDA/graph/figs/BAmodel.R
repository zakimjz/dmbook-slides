library(igraph)


cdeg = function(g){
	dV = sapply(V(g), function(v){ degree(g,v);});
	cdf = cumsum(dV)/sum(dV);
	return (cdf);
}


BAmodel = function(n0,q,t){
	#initial SW ring with k=1
	g <- watts.strogatz.game(1,n0,1,0);
	for (i in 1:t){
		cdf = cdeg(g);
		u = n0+i-1;
		g = add.vertices(g,1);
		e = 0;
		while (e < q){
			pv = runif(1,0,1);
			idx = which(cdf >= pv);
			v = head(idx,1);
			if (cdf[v] != pv){
				v= v-1; #if cdf matches, no need to decrement id
			}
			if (u != v && !are.connected(g,u,v)){
				g = add.edges(g,c(u,v));
				e = e+1;
			}
		}
	}
	return (g);
}

q=2
t=12
n0=3
g = BAmodel(n0,q,t)
l <- layout.circle(g)
R = c(0,-1,1,0) #rotate by -90deg
dim(R) = c(2,2)
l = l%*%R;
plot(g, layout=l)
l
(l+1)*2.5
E(g)


niters = 10;
n0=3;
q=3;
t=997;
maxdeg = 0;
dd = c(rep(0,n0+t+1));
dia = 0;
cc = 0;
for (iter in 1:niters){
	g = BAmodel(n0,q,t);
	dia = dia + diameter(g);
    cc = cc+ transitivity(g);
	degdist = degree.distribution(g);
	md = length(degdist);
	if (md > maxdeg){maxdeg = md;}
	for (i in 1:md){
		dd[i] = dd[i] + degdist[i];
	}
	cat("iter", iter, degdist, "\n");
}
dia/niters
cc/niters
dd[0:maxdeg]/niters

ndd = dd[0:maxdeg]/niters
degs = which(ndd>0)-1
ldd = sapply(ndd[ndd>0],function(x){log(x,2)})
ldegs = sapply(degs,function(x){log(x,2)})
ll = cbind(ldegs,ldd)
ll
lm(ldd[1:60]~ldegs[1:60])

