set xlabel "Radius (cm)"
set ylabel "Mass flow rate(g s^{-1})"

set logscale xy

plot for [i=0:100:2] "MdotvsR.txt" every :1::i::i