
set xlabel "Radius (cm)"
set ylabel "Time (days)"
set zlabel "Surface Density (g cm^{-2}" rotate parallel
#set hidden3d
unset logscale xyz
set grid
set logscale xz

splot "EvsR.txt" using 2:1:3 with lines