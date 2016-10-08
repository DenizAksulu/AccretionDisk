set xlabel "Time (day)"
set ylabel "Luminosity (erg/s)"


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

f(x) = j*x**8 + a*x**7 + b*x**6 + c*x**5 + d*x**4 + e*x**3 + f*x**2 + g*x + h
fit f(x) "lightcurve.txt" via a, b, c, d, e, f, g, h, j

plot f(x) with lines, "lightcurve.txt" with lines