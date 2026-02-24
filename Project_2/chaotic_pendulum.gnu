set terminal png size 800,600
set output 'poincare_section.png'
set title 'Poincaré Section - Driven Damped Pendulum'
set xlabel 'Theta (radians)'
set ylabel 'Omega (rad/s)'
set grid
plot 'poincare_section.out' using 1:2 with points pt 7 ps 0.1 lc rgb 'blue' notitle

set output 'phase_space_full.png'
set title 'Full Phase Space - Driven Damped Pendulum'
plot 'pendulum_trajectory.out' using 2:3 with points pt 7 ps 0.001 lc rgb 'red' title 'Full Trajectory'

set output 'time_series.png'
set title 'Time Series - Driven Damped Pendulum'
set xlabel 'Time (s)'
set ylabel 'Angle/Velocity'
plot 'pendulum_trajectory.out' using 1:2 with lines lw 1 title 'Theta (rad)', \
     'pendulum_trajectory.out' using 1:3 with lines lw 1 title 'Omega (rad/s)'
