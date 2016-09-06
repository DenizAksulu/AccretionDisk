set xlabel "Time (s)"
set ylabel "Luminosity (erg/s)"
samples(x) = $0 > 4 ? 5 : ($0+1)
avg5(x) = (shift5(x), (back1+back2+back3+back4+back5)/samples($0))
shift5(x) = (back5 = back4, back4 = back3, back3 = back2, back2 = back1, back1 = x)

#
# Initialize a running sum
#
init(x) = (back1 = back2 = back3 = back4 = back5 = sum = 0)
datafile = 'lightcurve.txt'
set logscale xy

plot "lightcurve.txt" with lines,\
	 sum = init(0), \
     datafile using 0:2 title 'data' lw 2 lc rgb 'forest-green', \
     '' using 0:(avg5($2)) title "running mean over previous 5 points" pt 7 ps 0.5 lw 1 lc rgb "blue", \
     '' using 0:(sum = sum + $2, sum/($0+1)) title "cumulative mean" pt 1 lw 1 lc rgb "dark-red"