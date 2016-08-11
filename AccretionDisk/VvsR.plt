set xlabel "Radius (cm)"
set ylabel "Viscosity"

set logscale xy

plot for [i=0:100:2] "VvsR.txt" every :1::i::i