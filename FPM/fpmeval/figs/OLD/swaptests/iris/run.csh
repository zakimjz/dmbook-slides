 perl swaptdb.pl iris.db irisswap 10 10
 perl swapfreq.pl iris.db irisswap irisswaptmp 10 10 10
 perl freqchanges.pl iris.db  iristmp1 iristmp2 10 10 10 10 > iris.out

awk '{print $1-$2, $0}' iris.out | sort -n -k1,1
awk '{print $1-$2, $0}' iris.out | sort -rn -k1,1
awk '{print $9-$10, $0}' iris.out | sort -n -k1,1
awk '{print $9-$10, $0}' iris.out | sort -rn -k1,1
awk '{print ($9-$10)/($9+0.000001), $0}' iris.out | sort -n -k1,1
awk '{print ($9-$10)/($9+0.000001), $0}' iris.out | sort -rn -k1,1

ls *.dat | awk '{print "sed -f cleanswap.sed " $1 "| awk -f statsfreq.awk > " $1".stat"}' | sh
cat irisswap.[0-9].dat.stat | awk -f sumstat.awk

