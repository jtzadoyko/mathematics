#!/usr/local/bin/gnuplot -persist
# set terminal x11  noraise nopersist
# set output

set key title ""
set key outside
set title "" 
set xlabel "" 
set ylabel "" 
set zero 1e-08

# Use together with line style "ls"
load "set1.pal"

plot [][] "data1.txt" u 1:10 w l ls 1

#    EOF
