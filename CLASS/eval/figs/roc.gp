set term postscript eps 20 enhanced;
set pointsize 2;
set xrange [0:1];
set yrange [0:1];

set key outside;
set grid;
set size square;
set xtics 0.2;
set ytics 0.2;
set xlabel 'FPR';
set ylabel 'TPR';

set output 'roc.eps';

plot 'roc.dat' u 1:2 w lp lw 5 title ''
