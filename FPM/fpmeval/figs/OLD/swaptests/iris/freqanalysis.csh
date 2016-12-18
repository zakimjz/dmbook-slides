./fim_all iris.db 1 iris.freq.1 > iris.freq.1.out
./fim_all iris.db 2 iris.freq.2 > iris.freq.2.out
./fim_all iris.db 3 iris.freq.3 > iris.freq.3.out
./fim_all iris.db 4 iris.freq.4 > iris.freq.4.out
./fim_all iris.db 5 iris.freq.5 > iris.freq.5.out
./fim_all iris.db 10 iris.freq.10 > iris.freq.10.out
./fim_all iris.db 15 iris.freq.15 > iris.freq.15.out
./fim_all iris.db 20 iris.freq.20 > iris.freq.20.out
./fim_all iris.db 25 iris.freq.25 > iris.freq.25.out
./fim_all iris.db 30 iris.freq.30 > iris.freq.30.out

sed 's/[()]//g' iris.freq.1 | awk -f statsfreq.awk  > iris.statsfreq.1
sed 's/[()]//g' iris.freq.2 | awk -f statsfreq.awk  > iris.statsfreq.2
sed 's/[()]//g' iris.freq.3 | awk -f statsfreq.awk  > iris.statsfreq.3
sed 's/[()]//g' iris.freq.4 | awk -f statsfreq.awk  > iris.statsfreq.4
sed 's/[()]//g' iris.freq.5 | awk -f statsfreq.awk  > iris.statsfreq.5
sed 's/[()]//g' iris.freq.10 | awk -f statsfreq.awk  > iris.statsfreq.10
sed 's/[()]//g' iris.freq.15 | awk -f statsfreq.awk  > iris.statsfreq.15
sed 's/[()]//g' iris.freq.20 | awk -f statsfreq.awk  > iris.statsfreq.20
sed 's/[()]//g' iris.freq.25 | awk -f statsfreq.awk  > iris.statsfreq.25
sed 's/[()]//g' iris.freq.30 | awk -f statsfreq.awk  > iris.statsfreq.30

sed 's/.*(//' iris.freq.1 | sed 's/)//' | sort -n | awk -f count.awk | sort -n -r | awk '{s+=$2; print $0, s/321.0}' > iris.freq.1.cdf
