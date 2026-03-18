set terminal gif animate delay 15 size 900,900
set output 'animation.gif'
set xrange [0:50]
set yrange [0:50]
set size square
set style fill solid 1.0
unset key
do for [i=0:749] {
    plot 'animation.dat' index i with points pt 7 ps 1.5 lc rgb 'blue'
}
