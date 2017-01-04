set xlabel "Radius (cm)"
set ylabel "Time (days)"
set zlabel "Irradiation Temperature(K)" rotate parallel
#set hidden3d
set logscale xyz
set grid
unset logscale yx
#set pm3d at st
#set logscale cb
splot "TirrvsR.txt" using 2:1:3 with lines