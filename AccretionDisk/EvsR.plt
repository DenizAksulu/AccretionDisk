set xlabel "Radius (cm)"
set ylabel "Surface Density (g cm^{-2})"

set logscale xy

plot for [i=0:100:20] "EvsR.txt" using 2:3 with lines, "EvsR_analytic.txt" using 2:3 with lines every:1::i::i 