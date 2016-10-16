set xlabel "Time (day)"
set ylabel "Hot radius (cm)"


samples(x) = $0 > 4 ? 5 : ($0+1)
avg5(y) = (shift5(y), (back1+back2+back3+back4+back5)/samples($0))
shift5(y) = (back5 = back4, back4 = back3, back3 = back2, back2 = back1, back1 = y)

#
# Initialize a running sum
#
init(x) = (back1 = back2 = back3 = back4 = back5 = sum = 0)
datafile = 'lightcurve.txt'
set logscale y
unset logscale x

plot "Rhot_T.txt" with lines