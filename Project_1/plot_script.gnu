set title 'Ping Pong Ball Trajectory'
set xlabel 'X (meters)'
set ylabel 'Y (meters)'
set zlabel 'Z (meters)'
set grid
splot 'pingpong_trajectory.out' with lines linewidth 2 title 'Trajectory'
pause -1 'Press any key to close'
