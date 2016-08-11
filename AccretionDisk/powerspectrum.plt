set xlabel "Energy (eV)"
set ylabel "Luminosity (erg s^{-1})"

set logscale xy

plot for [i=0:100:2] "powerspectrum.txt" every :1::i::i