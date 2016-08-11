set xlabel "Radius (cm)"
set ylabel "Surface Temperature(K)"

set logscale xy

plot for [i=0:100:2] "TsurvsR.txt" every :1::i::i