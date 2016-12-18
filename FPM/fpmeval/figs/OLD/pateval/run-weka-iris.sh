java -cp "./weka.jar"  weka.filters.supervised.attribute.Discretize -i iris.arff  -c last -o iris-discretizedMDL.arff
java -cp "./weka.jar" weka.associations.Apriori -I -N 100000 -T 1 -C 0.1 -M 0.065 -t iris-discretizedMDL.arff > iris-pats-rules-conf0.1-supp10.out
#find class rules
grep ": cl" iris.disc.weka.rules.fmt > iris.disc.weka.rules.cl
#sort by supp and conf
grep cl1 iris.disc.weka.rules.cl | sort -k3,5 -n | tail -1
grep cl2 iris.disc.weka.rules.cl | sort -k3,5 -n | tail -1
grep cl3 iris.disc.weka.rules.cl | sort -k3,5 -n | tail -1
#sort by conv and lift
grep cl1 iris.disc.weka.rules.cl | sort -k12,12 -k7,7 -n | tail -1
grep cl2 iris.disc.weka.rules.cl | sort -k12,12 -k7,7 -n | tail -1
grep cl3 iris.disc.weka.rules.cl | sort -k12,12 -k7,7 -n | tail -1
#rules containing pl2, pw2 and cl2 only
grep cl2 iris.disc.weka.rules.fmt | grep pl2 | grep pw2 | grep -v sl | grep -v sw
#print rule as points
awk '($0~"cl1"){print $3, $5}' iris.disc.weka.rules.cl > irissuppconfclass1.dat
awk '($0~"cl2"){print $3, $5}' iris.disc.weka.rules.cl > irissuppconfclass2.dat
awk '($0~"cl3"){print $3, $5}' iris.disc.weka.rules.cl > irissuppconfclass3.dat
awk '($0~"cl1"){print $7, $12}' iris.disc.weka.rules.cl > irisliftconvclass1.dat
awk '($0~"cl2"){print $7, $12}' iris.disc.weka.rules.cl > irisliftconvclass2.dat
awk '($0~"cl3"){print $7, $12}' iris.disc.weka.rules.cl > irisliftconvclass3.dat
#now do the pattern assessment via rule-based approach (avg lift)
#find all freq itemsets and rules with minsup=1
java -cp "./weka.jar" weka.associations.Apriori -I -N 100000 -T 1 -C 0.1 -M 0.005 -t iris-discretizedMDL.arff > iris-pats-rules-conf0.1-supp1.out
#there r 321 freq itemsets: F1(15), F2(72), F3(121), F4(89), F5(24)
#out of these the only itemsets with size >= 2 can contribute to rules,
#i.e., 306 patterns
awk '($7>=15&&$1>2.5&&$0~"cl")' result.case > result.case.supp15.lift2.5.cl
sort -k5,5 -n result.case.supp15.lift2.5.cl
