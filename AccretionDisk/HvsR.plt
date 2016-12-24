set xlabel "Radius (cm)"
set ylabel "Time (days)"
set zlabel "Scale Height (cm)" rotate parallel
#set hidden3d
unset logscale xyz
set grid
unset logscale y
#set pm3d at st
#set logscale cb
splot "HvsR.txt" using 2:1:3 with lines