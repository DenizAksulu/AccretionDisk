set xlabel "Radius (cm)"
set ylabel "Time (days)"
set zlabel "Shadow (1 if there is a shadow)" rotate parallel
#set hidden3d
unset logscale xyz
set grid
unset logscale y
#set pm3d at st
#set logscale cb
splot "ShadowsvsR.txt" using 2:1:3 with lines