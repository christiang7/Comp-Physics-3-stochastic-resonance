echo '----------------Start Programm----------------'

#### comiplation and producing of the output

clang++ sourcecode_output_Data/Main.cpp -o sourcecode_output_Data/stoch_reso && time ./sourcecode_output_Data/stoch_reso

#### Plot der Daten

echo '----------------Gnuplot Start----------------'
### plotting of the x values
gnuplot plotting.plt -p &

### plotting for task 2 plot kramers rule ln(T_K(D)) - 1/D
#gnuplot plotting-T_K-sigma.plt -p &


### plottting task 4
## plotting amplitude A_out over noise intensity D
#gnuplot plotting-ampl.plt -p &
## plotting maxima of A output over A input 
#gnuplot plotting-AO-AI.plt -p &

#gnuplot plotting-A_k-k.plt -p &
#gnuplot documentation/gnuplot_pictures/plot_ampl_phase.plt -p &
echo '----------------Ende ----------------'

