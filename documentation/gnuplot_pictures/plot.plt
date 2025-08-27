set term pdf

#    sigma = D
#set multiplot layout 2,2 title 'data for Amp=0.1 and omega=0.001 (<=> T=1000)'


  set xlabel "time t in s"
  set ylabel "x (a.u.)"

#set xtics 0, 1000, 3000
#set ytics -2.5, 1, 2.5


set xrange [0:2000]

# plot #1-4 here
set output 'Data_plot_D0.pdf'		
plot 'generated_data_for_plot/data_sigma0.000000.txt' pt 5 ps 0.03 title 'D=0'
set output 'Data_plot_D0p05.pdf'
plot 'generated_data_for_plot/data_sigma0.050000.txt' pt 5 ps 0.03 title 'D=0.05'
set output 'Data_plot_D0p1.pdf'
plot 'generated_data_for_plot/data_sigma0.100000.txt' pt 5 ps 0.03 title 'D=0.1'
set output 'Data_plot_D0p2.pdf'
plot 'generated_data_for_plot/data_sigma0.200000.txt' pt 5 ps 0.03 title 'D=0.20'



unset multiplot
exit
reset
