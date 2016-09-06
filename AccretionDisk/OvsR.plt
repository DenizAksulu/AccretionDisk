set xlabel "Radius (cm)"
set ylabel "Opacity"

set logscale xy

plot for [i=0:100:2] "OvsR.txt" every :1::i::i