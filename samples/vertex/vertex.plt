reset
set key
set zeroaxis lw 2

set style line 1 lt 1 lw 2 pt 6 ps 1.5
set style line 2 lt 2 lw 2 pt 4 ps 1.5
set style line 3 lt 3 lw 2 pt 6 ps 1.5
set style line 4 lt 4 lw 2 pt 4 ps 1.5

set style line 11 lt 1 lw 2 pt 1 ps 0.8 lc rgb "black"
set style line 12 lt 1 lw 2 pt 2 ps 0.8 lc rgb "black"
set style line 13 lt 1 lw 2 pt 3 ps 0.8 lc rgb "black"


set terminal postscript eps color enhanced "Times-Roman" 24
set output "vertex.eps"
set xlabel "{/Symbol e}_{/Times-Italic n}"
set ylabel "Re {/Symbol g}^{01}(i{/Symbol e}_{/Times-Italic n}, i{/Symbol e}_{{/Times-Italic n}{/Symbol \242}}, i{/Symbol w}_{/Times-Italic m})"
set xrange[*:*]
set yrange[*:*]

set key left bottom
set label "{/Symbol e}_{{/Times-Italic n}{/Symbol \242}} = {/Symbol p}{/Times-Italic T}" at graph 0.6, 0.05

plot "<awk '$2 == 4 && $3 == 0' vertex_lo.dat"u 4:9 title"(lo) {/Symbol w}_{/Times-Italic m} = 0" w lp ls 1\
, "<awk '$2 == 4 && $3 == 0' vertex_lo_G4.dat"u 4:9 title"" w p ls 11\
, "<awk '$2 == 4 && $3 == 0' vertex_tr.dat"u 4:9 title"(tr) {/Symbol w}_{/Times-Italic m} = 0" w lp ls 2\
, "<awk '$2 == 4 && $3 == 0' vertex_tr_G4.dat"u 4:9 title"" w p ls 12\

set output
