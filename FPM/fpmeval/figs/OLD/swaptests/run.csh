 perl swaptdb.pl simple.db simpleswap 1000 4
 perl swapfreq.pl simple.db simpleswap simpleswaptmp 2 1000 4
 perl freqchanges.pl simple.db  simpletmp1 simpletmp2 2 2 2 1000 > simple.out

 perl swaptdb.pl iris.db irisswap 10 10
 perl swapfreq.pl iris.db irisswap irisswaptmp 10 10 10
 perl freqchanges.pl iris.db  iristmp1 iristmp2 10 10 10 10 > iris.out
