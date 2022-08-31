reset
set key
set zeroaxis lw 2

set style line 1 lt 1 lw 2 pt 1 ps 1.5
set style line 2 lt 2 lw 2 pt 2 ps 1.5
set style line 3 lt 3 lw 2 pt 6 ps 1.5
set style line 4 lt 4 lw 2 pt 4 ps 1.5


set terminal postscript eps color enhanced "Times-Roman" 24
set output "Gf_t.eps"
set xlabel "{/Symbol t}"
set ylabel "{/Times-Italic G}_{imp}({/Symbol t})"
set xrange[0:*]
set yrange[:0]

plot "Gf_t.dat"u 1:2:3 title"0" w ye ls 1\
, ""u 1:4:5 title"1" w ye ls 2\

set output


set terminal postscript eps color enhanced "Times-Roman" 24
set output "Gf_w.eps"
set xlabel "{/Symbol e}_{/Times-Italic n}"
set ylabel "{/Times-Italic G}_{imp}(i{/Symbol e}_{/Times-Italic n})"
set xrange[0:10]
set yrange[*:*]

set key bottom

plot "Gf_w.dat"u 1:2 title"0 Re" ls 1, ""u 1:3 title"0 Im" ls 2\
, ""u 1:4 title"1 Re" ls 3, ""u 1:5 title"1 Im" ls 4\

set output


set terminal postscript eps color enhanced "Times-Roman" 24
set output "Gf_pade.eps"
set xlabel "{/Symbol w}"
set ylabel "{/Symbol r}_{imp}({/Symbol w})"
set xrange[*:*]
set yrange[0:*]

set key top

plot "Gf_pade.dat"u 1:(-$3/pi) title"0" w l ls 1\
, ""u 1:(-$5/pi) title"1" w l ls 2

set output



set terminal postscript eps color enhanced "Times-Roman" 24
set output "self_w.eps"
set xlabel "{/Symbol e}_{/Times-Italic n}"
set ylabel "{/Symbol S}_{imp}(i{/Symbol e}_{/Times-Italic n})"
set xrange[0:10]
set yrange[*:*]

set key bottom

plot "self_w.dat"u 1:2 title"0 Re" ls 1, ""u 1:3 title"0 Im" ls 2\
, ""u 1:4 title"1 Re" ls 3, ""u 1:5 title"1 Im" ls 4\

set output



set terminal postscript eps color enhanced "Times-Roman" 24
set output "self_w_dyson.eps"
set xlabel "{/Symbol e}_{/Times-Italic n}"
set ylabel "{/Symbol S}_{imp}(i{/Symbol e}_{/Times-Italic n})"
set xrange[0:10]
set yrange[*:*]

set key bottom

plot "self_w_dyson.dat"u 1:2 title"0 Re" ls 1, ""u 1:3 title"0 Im" ls 2\
, ""u 1:4 title"1 Re" ls 3, ""u 1:5 title"1 Im" ls 4\

set output



set terminal postscript eps color enhanced "Times-Roman" 24
set output "GSigma_t.eps"
set xlabel "{/Symbol t}"
set ylabel "({/Times-Italic G}{/Symbol S})({/Symbol t})"
set xrange[0:*]
set yrange[:0]

plot "GSigma_t.dat"u 1:2:3 title"0" w ye ls 1\
, ""u 1:4:5 title"1" w ye ls 2\

set output
