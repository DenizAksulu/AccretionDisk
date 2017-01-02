set xlabel "Time (day)"
set ylabel "Luminosity (erg/s)"

set logscale y
unset logscale x

plot[:][1e10:] "lightcurve_optical.txt" with lines, "lightcurve_HEXTE.txt" with lines