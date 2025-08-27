reset
#set key outside rmargin
set grid
set xrange [0:12]
set ylabel 'A_{max}'
set term pdf color enhanced
###
set xlabel 'A_{in}'
set output sprintf('plotting-AO-AI.pdf')
f(x)=a/(x+b)
fit f(x) 'maxima-A-out.txt' using ($1):(($3)/($1)) via a,b
plot  f(x) with lines lw 1.5 title 'A_{fit}(A_{in})', 'maxima-A-out.txt' using ($1):(($3)/($1)) title 'A_{max}(A_{in})' with points pt 5 ps 0.5

set term qt
replot

