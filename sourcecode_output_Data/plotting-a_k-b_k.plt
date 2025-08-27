reset
set grid
#set yrange [-2:3]
set ylabel 'T_K'
#set term png
#set logscale y
###
#set xrange [0:100]
set xlabel '1/sigma/sigma'
#set arrow from 5.70477,0 to 5.70477,2 nohead lw 2
#set output sprintf('plotting-A-h0.png')
T_K(x) = log(a*x+ b) 
fit T_K(x) 'T_K.txt' using ($1):($2) via a,b
plot "T_K.txt" using (($1)):(($2)) lt rgb 'red' title 'T_K', T_K(x) title "Fit Function" #

#set term qt
#replot
