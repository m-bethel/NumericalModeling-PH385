set title 'Mean Squared Displacement vs Step'
set xlabel 'Step'
set ylabel '<r^2>'
set datafile separator ","
set key noautotitle
plot 'msd.csv' skip 1 using 1:2 with lines linewidth 2
pause -1 'Press any key to close'
