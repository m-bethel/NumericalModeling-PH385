set term qt 1 title 'Temperature vs E' noraise
set datafile separator ','
set terminal png size 1200,600
set output 'phase_transition.png'
set xlabel 'Total Internal Energy (E)'
set ylabel 'Temperature (k_B T)'
set grid
set key left top
set offsets graph 0.05, 0.05, 0.05, 0.05
set title 'Phase Transition'
plot 'heating.csv' using 1:2 with lines lc rgb 'red' lw 1 title 'Heating', \
     'cooling.csv' using 1:2 with lines lc rgb 'blue' lw 1 title 'Cooling'
