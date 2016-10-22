set xlabel "Radius (cm)"
set ylabel "Shadow (1 if there is a shadow)"

unset logscale y
set logscale x

plot for [i=0:100:2] "ShadowvsR.txt" every :1::i::i with points pointsize 10