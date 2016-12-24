set xlabel "Radius (cm)"
set ylabel "Time (days)"
set zlabel "Surface Density (g cm^{-2})" rotate parallel
#set hidden3d
set logscale xyz
set grid
unset logscale y
#set pm3d at st
#set logscale cb
splot "EvsR.txt" using 2:1:3 with lines