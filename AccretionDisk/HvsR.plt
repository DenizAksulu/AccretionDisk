set xlabel "Radius (cm)"
set ylabel "Scale Height (cm)"
set yrange [0.001:1e10]
set logscale xy

plot for [i=0:100:2] "HvsR.txt" every :1::i::i