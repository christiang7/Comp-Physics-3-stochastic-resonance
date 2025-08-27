reset
set key outside rmargin
set grid
#set yrange [0:2]
set ylabel 'A_{out}/A_{in}'
set term pdf color enhanced
###
#set xrange [0:2.2]
set xlabel 'D'
#set arrow from 5.70477,0 to 5.70477,2 nohead lw 2
set output sprintf('plotting-A-D-all.pdf')
plot for [i=3:5] 'documentation/gnuplot_pictures/A_k-phi_k_omega0.001000_A0.0'.i.'.txt' using ($1):($2/(i*0.01)) title 'A_{in} = 0.0'.i with lines, for [i=10:60:5] 'documentation/gnuplot_pictures/A_k-phi_k_omega0.001000_A0.'.i.'.txt' using ($1):($2/(i*0.01)) title 'A_{in} = 0.'.i with lines, for [i=1:4] 'documentation/gnuplot_pictures/A_k-phi_k_omega0.001000_A'.i.'.00.txt' using ($1):($2/(i)) title 'A_{in} = '.i.'.0' with lines, 'documentation/gnuplot_pictures/A_k-phi_k_omega0.001000_A10.00.txt' using ($1):($2/(10)) title 'A_{in} = 10.0' with lines

#plot "b_x_ampl-pl.txt" index 0 using 1:2 title 'x' lt rgb 'blue' with lines, "b_x_ampl-pl.txt" index 1 using 1:2 title 'x' lt rgb 'blue' with lines, "b_x_ampl-pl.txt" index 2 using 1:2 title 'x' lt rgb 'blue' with lines, "b_x_ampl-pl.txt" index 3 using 1:2 title 'x' lt rgb 'blue' with lines, "b_x_ampl-pl.txt" index 4 using 1:2 title 'x' lt rgb 'blue' with lines, "b_x_ampl-pl.txt" index 5 using 1:2 title 'x' lt rgb 'blue' with lines

set term qt
replot

