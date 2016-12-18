#!/bin/bash
#do S swaps, generate D datasets
S=150
D=100
perl swaptdb.pl simple.db simple/simpleswap $S $D
rm simple.tot.out
for ((i=0; i<=$D; i++))
do
  #mine with minsup 3
  echo "process: $i"
  ./fim_all simple/simpleswap.$i.dat 3 simple/simpleswap.$i.freq
  sort -n simple/simpleswap.$i.freq -o simple/simpleswap.$i.freq
  sed 's/[()]//g' simple/simpleswap.$i.freq | awk -f statsfreq.awk > simple/simpleswap.freq.$i.stats
  tail -5 simple/simpleswap.freq.$i.stats
  tail -1 simple/simpleswap.freq.$i.stats >> simple.tot.out
  diff simple/simpleswap.0.freq simple/simpleswap.$i.freq
done
#check out stats for 1245
echo "stats of 1245"
#awk '{if (NF>4) print $0}' simple/simpleswap.*.freq | grep 1 | grep 2 | grep 4| grep 5 > simple.1245.out
awk '{if (NF>4) print $0}' simple/simpleswap.*.freq | perl -wne 'chop; print sort split; print "\n"' | grep "1245" > simple.1245.out
wc simple.1245.out
echo "stats of 234"
#awk '{if (NF==4) print $0}' simple/simpleswap.*.freq | grep 2 | grep 3| grep 4 > simple.234.out
awk '{if (NF==4) print $0}' simple/simpleswap.*.freq | perl -wne 'chop; print sort split; print "\n"' | grep "234" > simple.234.out
wc simple.234.out
