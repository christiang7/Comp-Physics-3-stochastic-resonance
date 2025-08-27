reset
set grid
set ylabel 'A'
set boxwidth 0.01
set style fill solid
plot "A_k-phi_k.txt" using 1:2 title 'A_k' with boxes

