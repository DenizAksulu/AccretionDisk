set xlabel "Radius (cm)"
set ylabel "Pressure"

set logscale xy

plot for [i=0:100:2] "PgasvsR.txt" every :1::i::i