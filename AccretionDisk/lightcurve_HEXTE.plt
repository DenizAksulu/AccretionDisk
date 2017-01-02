set xlabel "Time (day)"
set ylabel "Luminosity (erg/s)"

set logscale y
unset logscale x

plot[:][1e5:] "lightcurve_HEXTE.txt" with lines