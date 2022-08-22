reset
set key
set zeroaxis lw 2

set style line 1 lt 1 lw 2 pt 1 ps 1.5
set style line 2 lt 2 lw 2 pt 2 ps 1.5
set style line 3 lt 3 lw 2 pt 6 ps 1.5
set style line 4 lt 4 lw 2 pt 4 ps 1.5


set terminal postscript eps color enhanced "Times-Roman" 24
set output "delta_tau.eps"
set xlabel "{/Symbol t}"
set ylabel "{/Symbol D}({/Symbol t})"
set xrange[0:*]
set yrange[:0]

plot "delta_tau.dat"u 1:2 title"0" w lp ls 1\
, ""u 1:3 title"1" w lp ls 2

set output



set terminal postscript eps color enhanced "Times-Roman" 24
set output "delta_omega.eps"
set xlabel "{/Times-Italic n}"
set ylabel "{/Symbol D}(i{/Symbol e}_{/Times-Italic n})"
set xrange[0:10]
set yrange[*:*]

set key bottom

plot "delta_omega.dat"u 1:2 title"0 real" ls 1, ""u 1:3 title"0 imag" ls 2\
#, ""u 1:4 title"1 real" ls 3, ""u 1:5 title"1 imag" ls 4\

set output
