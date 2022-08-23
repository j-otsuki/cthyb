reset
set key


set style line 1 lt 1 lw 2 pt 1 ps 1.5
set style line 2 lt 2 lw 2 pt 2 ps 1.5
set style line 3 lt 3 lw 2 pt 6 ps 1.5
set style line 4 lt 4 lw 2 pt 4 ps 1.5


set terminal postscript eps color enhanced "Times-Roman" 24
set output "stat.eps"
set xlabel "{/Times-Italic k}"
set ylabel ""
set xrange[0:*]
set yrange[0:]

plot "stat.dat"u 1:2 title"0" w lp ls 1\
, ""u 1:3 title"1" w lp ls 2\
#, ""u 1:4 title"2" w lp ls 3\

set output
