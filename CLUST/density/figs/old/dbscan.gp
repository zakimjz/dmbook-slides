set term postscript eps enhanced "Times-BoldItalics" 24
set output 'dbscan.eps'

set noxtics
set noytics

plot 'dbscan.dat' u 1:2 w p pt 6 ps 1.5 t''
