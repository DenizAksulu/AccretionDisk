set xlabel "Radius (cm)"
set ylabel "Central Temperature(K)"

set logscale xy

plot for [i=0:100:2] "TcvsR.txt" every :1::i::i