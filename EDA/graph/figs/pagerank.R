norm = function(x){
	sqrt(sum(x*x));
}

myeig = function(X,errt){
	p = c(1,1,1,1,1)
	oldeig = 1;
	iter = 0;
	while(T){
        np = X %*% p;
        mm = max(np);
        ii = which(np == mm);
        neweig = np[ii[1]]/p[ii[1]];
        cat(iter, " ", ii[1], " ", neweig,"\n");
		cat(np, " ", np/mm, "\n");
        err = abs(neweig-oldeig);
		diff = np/mm-p;
        err = sum(abs(diff))
		
        if (err < errt){
                cat("err", err, "\n");
                break;
        }
        p = np/mm;
        iter = iter+1;
        oldeig=neweig;
	}
	p = p/(sqrt(sum(p*p)))
}

myeigUNnormalized = function(X,iters){
	p = c(1,1,1,1,1)
	oldeig = 1;
	iter = 0;
	while(T){
        np = X %*% p;
        mm = max(np);
        ii = which(np == mm);
        neweig = np[ii[1]]/p[ii[1]];
        cat(iter, " ", ii[1], " ", neweig,"\n");
		cat(np, "\n");
        if (iter > iters){
            break;
        }
        p = np;
        iter = iter+1;
        oldeig=neweig;
	}
	p = p/(sqrt(sum(p*p)))
}

hits <- function(X, errt){
	a = c(1,1,1,1,1)
	h = c(1,1,1,1,1)
	oldeig = 1;
	iter = 0;
	while(T){
        nh = X %*% a;
		mh = max(nh);
        ih = which(nh == mh);
		nh = nh/ih[1];
		na = t(X) %*% nh;
		#na = na/norm(na);
		ma = max(na);
        ia = which(na == ma);
		na = na/ia[1];
		diff = (na-a)+(nh-h);
        err = sum(abs(diff))
        if (err < errt){
                cat("err", err, "\n");
                break;
        }
        iter = iter+1;
		a = na;
		h = nh;
	}
	a = a/norm(a);
	h = h/norm(h);
	return (list(a=a,h=h));
}

setwd('/Users/zaki/research/DataMiningBook/dm08/EDA/graph/figs')
A = matrix(c(0,0,1,0,0,0,0,0,1,1,0,1,0,1,0,1,0,0,0,0,0,1,0,1,0),
		nrow=5, ncol=5)

print(A)
N =
matrix(c(0,0,1,0,0,0,0,0,1/3,1,0,0.5,0,1/3,0,1,0,0,0,0,0,0.5,0,1/3,0),
		nrow=5, ncol=5)

Nf=
matrix(c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2),
		nrow= 5, ncol=5)
d = 0.1
M = (1-d)*N + d*Nf
M
ee =eigen(t(M))
ee
cat("MYEIG-A\n")
p = myeig(t(A), 0.0001)
p
cat("MYEIGU-A\n")
p = myeigUNnormalized(t(A), 16)
p

cat("MYEIG-M\n")
p = myeig(t(M), 0.0001)
p

cat("HITS\n")
hh = hits(A,0.0001)
hh$a
hh$h
