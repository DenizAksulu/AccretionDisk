set xlabel "Radius (cm)"
set ylabel "Time (days)"
set zlabel "Alpha parameter" rotate parallel
#set hidden3d
unset logscale xyz
set grid
unset logscale y
#set pm3d at st
#unset logscale cb
splot "AvsR.txt" using 2:1:3 with lines #palette 