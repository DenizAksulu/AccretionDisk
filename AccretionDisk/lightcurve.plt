set xlabel "Time (day)"
set ylabel "Luminosity (erg/s)"

set logscale y
unset logscale x

plot[:][1e33:] "lightcurve.txt" with lines, "lightcurve_bolo.txt" with lines