#!/usr/bin/env bash
#CDF for frequent itemsets
./fim_all iris.db 1 iris.freq.1 > iris.freq.1.out
sed 's/.*(//' iris.freq.1 | sed 's/)//' | sort -n | awk -f iris/count.awk | sort -n | awk '{s+=$2; print $0, s/321.0}' > iris.freq.1.cdf
awk '{print "{"$1 ", " $3 "}"}' iris.freq.1.cdf > iris-supp-1.cdf
#
#now do relative lift analysis
S=2250 # steps 150 rows x 15 items
D=100 # number of swapped datasets
msup=10 #min sup
for ((i=0; i<=$D; i++))
do
perl freqchanges.pl iris.db iris/iristmp1 iris/iristmp2 $msup $msup $msup $S > iris/iris.$i.out
done
./rellift_iris.py iris/iris.*.out > iris.rellift.out
awk '{s+=1; print $0, s/140.0}' iris.rellift.out | awk '{print $1 ", " $NF }'
awk '{print ($1-$2)/$1, " ", ($9-$10)/$9 }' iris/iris.1.out
sed -f deconv.sed iris.rellift.out
