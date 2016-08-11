set xlabel "Radius (cm)"
set ylabel "Irradiation Temperature(K)"

set logscale xy

plot for [i=0:100:2] "TirrvsR.txt" every :1::i::i