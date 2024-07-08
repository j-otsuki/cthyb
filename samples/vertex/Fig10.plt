reset
set key
set zeroaxis lw 2

set style line 1 lt 1 lw 2 pt 7 ps 1.0 lc 1
set style line 2 lt 1 lw 2 pt 5 ps 1.0 lc 2
set style line 3 lt 1 lw 2 pt 9 ps 1.0 lc 3
set style line 11 lt 1 lw 2 pt 1 ps 0.8 lc rgb "black"
set style line 12 lt 1 lw 2 pt 2 ps 0.8 lc rgb "black"
set style line 13 lt 1 lw 2 pt 3 ps 0.8 lc rgb "black"




set terminal postscript eps color enhanced "Times-Roman" 24
set output "Fig10.eps"
set xlabel "{/Symbol e}_{/Times-Italic n}"
set ylabel "Re {/Symbol g}^{uu}(i{/Symbol e}_{/Times-Italic n}, i{/Symbol e}_{/Times-Italic n'}, i{/Symbol w}_{/Times-Italic m})"
set xrange[-20:20]
set yrange[-6:7]

set key left bottom

plot "<awk '$2 == 4 && $3 == 0' vertex2_lo.dat"u 4:(-$7) title"" w p ls 11\
, "<awk '$2 == 4 && $3 == 2' vertex2_lo.dat"u 4:(-$7) title"" w p ls 12\
, "<awk '$2 == 4 && $3 == 10' vertex2_lo.dat"u 4:(-$7) title"" w p ls 13\
, "<awk '$2 == 4 && $3 == 0' vertex_lo.dat"u 4:(-$7) title"{/Symbol w}_{/Times-Italic m} = 0" w lp ls 1\
, "<awk '$2 == 4 && $3 == 2' vertex_lo.dat"u 4:(-$7) title"4{/Symbol p}{/Times-Italic T}" w lp ls 2\
, "<awk '$2 == 4 && $3 == 10' vertex_lo.dat"u 4:(-$7) title"20{/Symbol p}{/Times-Italic T}" w lp ls 3\

set output
