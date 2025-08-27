set term pdf
set output 'ampl_1p0_phase--.pdf'


  set xlabel "D"
  set ylabel "Amplitude A (a.u.)"
  set y2label "Phase {/Symbol f}"

set ytics nomirror
set y2tics

set y2range [0:2]

path='A_k-phi_k_omega0.001000_A0.100000.txt'

		
#plot path u 1:2:4 with yerrorbars pt 5 ps 0.2 title 'Amplitude' axes x1y1 , path u 1:3:5 with yerrorbars pt 5 ps 0.2 title 'Phase' axes x1y2 

plot path u 1:2 pt 5 ps 0.2 title 'Amplitude' axes x1y1 , path u 1:3 pt 5 ps 0.2 title 'Phase' axes x1y2 
set term qt
replot


set term pdf
set output 'phase.pdf'
unset y2tics
unset y2label
plot path u 1:6 title "a_k Amplitude" pt 5 ps 0.2 lc 3, path u 1:7 title "b_k Amplitude" pt 5 ps 0.2 lc 7
set term qt
replot
exit
reset
