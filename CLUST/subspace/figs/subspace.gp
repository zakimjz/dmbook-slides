set term postscript eps enhanced "Times-Bold" 24
set xlabel "X"; 
set ylabel "Y";
set zlabel "Z";

set xrange [0:5]
set yrange [0:5]
set zrange [0:5]
set ticslevel 0;
set view 70,140
set border 120;
set style line 20 lt 1 lw 5
set grid xtics ytics ztics

set pointsize 2

set output 'subspace.eps';
set label 1 "C_1" at 0,0.25,4.4
set label 2 "C_2" at 0,4,2.4
splot 'subspace.dat' i 0 w l ls 20 t '', '' i 1 w l ls 20 t '',\
'' i 2 w l ls 20 t '',\
'' i 3 w l ls 20 t '', '' i 4 w l ls 20 t '',\
'' i 5 w l ls 20 t '', '' i 6 w l ls 20 t '',\
'' i 7 w l ls 20 t '', '' i 8 w l ls 20 t '',\
'' i 9 w l ls 20 t '', '' i 10 w l ls 20 t '',\
'' i 11 w l ls 20 t '',\
'' i 12 w p pt 6 t '', '' i 13 w p pt 6 t ''

set output 'subspace-xy.eps'
unset label 1 
unset label 2
set xlabel "X"
set ylabel "Y"
set label 1 "C_1" at 0.5,1.25
set label 2 "C_2" at 0.5,2.75
plot 'subspace.dat' i 12 u 1:2 w p pt 6 t '',\
'' i 13 u 1:2 w p pt 6 t '',\
'' i 0 u 1:2 w l ls 20 t '', '' i 1 u 1:2 w l ls 20 t '',\
'' i 2 u 1:2 w l ls 20 t '',\
'' i 3 u 1:2 w l ls 20 t '', '' i 4 u 1:2 w l ls 20 t '',\
'' i 5 u 1:2 w l ls 20 t '', '' i 6 u 1:2 w l ls 20 t '',\
'' i 7 u 1:2 w l ls 20 t '', '' i 8 u 1:2 w l ls 20 t '',\
'' i 9 u 1:2 w l ls 20 t '', '' i 10 u 1:2 w l ls 20 t '',\
'' i 11 u 1:2 w l ls 20 t ''

set output 'subspace-xz.eps'
unset label 1 
unset label 2
set xlabel "X"
set ylabel "Z"
set label 1 "C_1" at 1.25,3.1
set label 2 "C_2" at 4.25,1.1
plot 'subspace.dat' i 12 u 1:3 w p pt 6 t '',\
'' i 13 u 1:3 w p pt 6 t '',\
'' i 0 u 1:3 w l ls 20 t '', '' i 1 u 1:3 w l ls 20 t '',\
'' i 2 u 1:3 w l ls 20 t '',\
'' i 3 u 1:3 w l ls 20 t '', '' i 4 u 1:3 w l ls 20 t '',\
'' i 5 u 1:3 w l ls 20 t '', '' i 6 u 1:3 w l ls 20 t '',\
'' i 7 u 1:3 w l ls 20 t '', '' i 8 u 1:3 w l ls 20 t '',\
'' i 9 u 1:3 w l ls 20 t '', '' i 10 u 1:3 w l ls 20 t '',\
'' i 11 u 1:3 w l ls 20 t ''

set output 'subspace-yz.eps'
set xlabel "Y"
set ylabel "Z"
unset label 1 
unset label 2
set label 1 "C_1" at 0.25,4.4
set label 2 "C_2" at 4,2.4
plot 'subspace.dat' i 12 u 2:3 w p pt 6 t '', '' i 13 u 2:3 w p pt 6 t '',\
'' i 0 u 2:3 w l ls 20 t '', '' i 1 u 2:3 w l ls 20 t '',\
'' i 2 u 2:3 w l ls 20 t '',\
'' i 3 u 2:3 w l ls 20 t '', '' i 4 u 2:3 w l ls 20 t '',\
'' i 5 u 2:3 w l ls 20 t '', '' i 6 u 2:3 w l ls 20 t '',\
'' i 7 u 2:3 w l ls 20 t '', '' i 8 u 2:3 w l ls 20 t '',\
'' i 9 u 2:3 w l ls 20 t '', '' i 10 u 2:3 w l ls 20 t '',\
'' i 11 u 2:3 w l ls 20 t ''
