set term postscript eps 20 enhanced;
set pointsize 2;
set xrange [0:10];
set yrange [0:10];

set key outside;
set grid;
set size square;
set xtics 1;
set ytics 1;
set xlabel 'X_1';
set ylabel 'X_2';

set output 'svm-hull.eps';
plot 'svm.dat' index 0 w p pt 8 title '', \
     'svm.dat' index 1 w p pt 6 title '', \
      'svm.dat' index 4 w p pt 8 title '',\
      'svm.dat' index 5 w p pt 6 title '',\
      'svm.dat'	index 6 w l lt 3 title '',\
      'svm.dat' index 7 w l lt 3 title '',\
      4*x-22 title '' w l lt 2,\
      4*x-8 title '' w l lt 3,\
      4*x-15 title 'h' w l lt 1 lw 3;
