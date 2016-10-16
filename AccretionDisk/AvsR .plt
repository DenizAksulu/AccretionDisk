set xlabel "Radius (cm)"
set ylabel "Alpha parameter"
unset logscale xy

plot for [i=0:100:2] "AvsR.txt" every :1::i::i