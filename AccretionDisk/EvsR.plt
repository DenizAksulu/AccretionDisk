set xlabel "Radius (cm)"
set ylabel "Surface Density (g cm^{-2})"

set logscale y
unset logscale x

plot for [i=0:100:2] "EvsR.txt" every :1::i::i