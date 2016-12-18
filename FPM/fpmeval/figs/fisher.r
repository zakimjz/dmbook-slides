#!/usr/bin/env Rscript
library(stringr)
setwd('/Users/zaki/research/DataMiningBook/dm08/FPM/fpmeval/figs')

runcmd = function(cmd){
  res = system(cmd, intern=TRUE)
  res=str_trim(res)
  ca = strsplit(res, ' ')
  v = as.integer(unlist(ca)[1])
  return(v)
}

getP = function(aa,bb,cc,dd){
#  v = (factorial(aa+bb)*
#	  factorial(cc+dd)*
#	  factorial(aa+cc)*
#	  factorial(bb+dd))/
#	  (factorial(aa+bb+cc+dd)*
#	   factorial(aa)*
#	   factorial(bb)*
#	   factorial(cc)*
#	   factorial(dd))
  n = aa+bb+cc+dd
  v = choose(aa+bb, aa)*choose(cc+dd,cc)/choose(n, aa+cc)
  print(c(aa,bb,cc,dd,v))
  return(v)
}

getpval = function(aa,bb,cc,dd){
  pvalue = 0
  for (i in 0:min(bb,cc)){
	print(i)
	pvalue = pvalue + getP(aa+i, bb-i, cc-i, dd+i)
  }
  return (pvalue)
}

#######main
aa = runcmd("grep cl2 iris-discrete.txt | grep pw2 | wc")
bb = runcmd("grep -v cl2 iris-discrete.txt | grep pw2 | wc")
cc = runcmd("grep cl2 iris-discrete.txt | grep -v pw2 | wc")
dd = runcmd("grep -v cl2 iris-discrete.txt | grep -v pw2 | wc")
pvalue = getpval(aa,bb,cc,dd)
print(c("pw2 cl2", aa, bb, cc, dd, pvalue))

#now do the test for sw1, pw2, cl2
#R1: pw2 -> cl2
system("grep pw2 iris-discrete.txt > /tmp/pw2.txt")
aa = runcmd("grep cl2 /tmp/pw2.txt | grep sw1 | wc")
bb = runcmd("grep -v cl2 /tmp/pw2.txt | grep sw1 | wc")
cc = runcmd("grep cl2 /tmp/pw2.txt | grep -v sw1 | wc")
dd = runcmd("grep -v cl2 /tmp/pw2.txt | grep -v sw1 | wc")
pvalue = getpval(aa,bb,cc,dd)
print(c("sw1 cl2 | pw2", aa, bb, cc, dd, pvalue))

#R2: sw12 -> cl2
system("grep sw1 iris-discrete.txt > /tmp/sw1.txt")
aa = runcmd("grep cl2 /tmp/sw1.txt | grep pw2 | wc")
bb = runcmd("grep -v cl2 /tmp/sw1.txt | grep pw2 | wc")
cc = runcmd("grep cl2 /tmp/sw1.txt | grep -v pw2 | wc")
dd = runcmd("grep -v cl2 /tmp/sw1.txt | grep -v pw2 | wc")
pvalue = getpval(aa,bb,cc,dd)
print(c("pw2 cl2 | sw1", aa, bb, cc, dd, pvalue))

#R3: {} -> cl2
aa = runcmd("grep cl2 iris-discrete.txt | grep pw2 | grep sw1 | wc")
bb = runcmd("grep -v cl2 iris-discrete.txt | grep pw2 | grep sw1 | wc")
cc1 = runcmd("grep cl2 iris-discrete.txt | grep -v pw2 | wc")
cc2 = runcmd("grep cl2 iris-discrete.txt | grep -v sw1 | wc")
cc3 = runcmd("grep cl2 iris-discrete.txt | grep -v pw2 | grep -v sw1 | wc")
cc = cc1+cc2-cc3 #inclusion-exclusion
dd1 = runcmd("grep -v cl2 iris-discrete.txt | grep -v pw2 | wc")
dd2 = runcmd("grep -v cl2 iris-discrete.txt | grep -v sw1 | wc")
dd3 = runcmd("grep -v cl2 iris-discrete.txt | grep -v pw2 | grep -v sw1 | wc")
dd = dd1+dd2-dd3 #inclusion-exclusion
pvalue = getpval(aa,bb,cc,dd)
print(c("{} | sw1", aa, bb, cc, dd, pvalue))


