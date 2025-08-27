reset
set grid
#set yrange [-2:3]
set ylabel 'x'
#set term png
###
#set xrange [0:1.2]
set xlabel 't'
#set arrow from 5.70477,0 to 5.70477,2 nohead lw 2
#set output sprintf('plotting-A-h0.png')
plot "data.txt" using ($1):($2) title 'x' lt rgb 'blue' with lines #, "T_K0.txt" using ($1):($2) lt rgb 'green' title 'T_K0'#, "T_K.txt" using ($1):($2) lt rgb 'red' title 'T_K' #

#set term qt
#replot
