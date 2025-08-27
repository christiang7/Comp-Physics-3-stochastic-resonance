reset
set term pdf color enhanced
set key left
set grid
#set yrange [-2:3]
set ylabel 'ln(T_K)(D)'
#set logscale y
###
#set xrange [0:100]
set xlabel '1/D'
#set arrow from 5.70477,0 to 5.70477,2 nohead lw 2
#set output sprintf('plotting-A-h0.png')
#T_K(x) = log(a*x+ b) 
#T_K(x) = sqrt(a*x)+ b
T_K(x) = a*x+ b
f(x) = 0.25*x+log(sqrt(2)*pi)
fit T_K(x) 'T_K.txt' using ($1):(($2)) via a,b
set output sprintf('plotting-T_K.pdf')
plot "T_K.txt" using (($1)):(($2)) lt rgb 'red' pt 5 ps 0.2 title 'data ln(T_d)(D)',T_K(x) lw 1.5 title "Fit Function ln(T_f)(D)", f(x) lw 1.5 title 'Kramer formula ln(T_K)(D)'#

set term qt
replot
