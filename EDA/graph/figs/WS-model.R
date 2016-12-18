library(igraph)
g <- watts.strogatz.game(1, 8, 2, 0.0)
l <- layout.circle(g)
plot(g, layout=l)
l
(l+1)*2.5

g <- watts.strogatz.game(1, 20, 3, 0)
l <- layout.circle(g)
plot(g, layout=l)
l
(l+1)*2.5



edeg2 <- function(n,k,r){ 
	t2 = (3*(k-1))/((1+r)*(4*k*r+2*(2*k-1)));
	return (t2);
}


edeg <- function(n,k,r){ 
	t1 = (2*k*r/(n-2*k-1));
	t2 = (3*(k-1))/((1+r)*(4*k*r+2*(2*k-1)));
	return (t1+t2);
}

n = 1000;
k = 3;
runs = 40;
iters = 10;
X = matrix(0,(runs+1),6)
for (iter in 1:iters){
	for (i in 0:runs){
		r = i*0.005;
		g <- watts.strogatz.game(1, n, k, 0);
		mnew = as.integer(k*n*r);
		if (mnew > 0){
			for (nm in 1:mnew){
				e= sample(0:(n-1),2);
				g = add.edges(g,e);
			}
		}
		dd = diameter(g);
		cc = transitivity(g);
		pcc = edeg2(n,k,r);
		cat(r, mnew, dd, cc, pcc, "\n");
		X[(i+1),] = X[(i+1),] + c(r,mnew,dd,cc,pcc,pcc2);
	}
}
X = X/iters
X
