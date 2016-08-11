set xlabel "Radius (cm)"
set ylabel "Effecive Temperature(K)"

set logscale xy

plot for [i=0:100:2] "TeffvsR.txt" every :1::i::i