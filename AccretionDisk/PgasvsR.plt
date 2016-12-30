set xlabel "Radius (cm)"
set ylabel "Time (days)"
set zlabel "Pressure" rotate parallel
#set hidden3d
set logscale xyz
set grid
unset logscale y
#set pm3d at st
#set logscale cb
splot "PgasvsR.txt" using 2:1:3 with lines