set xlabel "Radius (cm)"
set ylabel "Gas/Radiation"

set logscale xy

plot for [i=1:1000:5] "PGvsR.txt" every :1::i::i