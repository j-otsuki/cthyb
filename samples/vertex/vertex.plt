reset
set key
set zeroaxis lw 2

set style line 1 lt 1 lw 2 pt 1 ps 1.5
set style line 2 lt 2 lw 2 pt 2 ps 1.5
set style line 3 lt 3 lw 2 pt 6 ps 1.5
set style line 4 lt 4 lw 2 pt 4 ps 1.5



set terminal postscript eps color enhanced "Times-Roman" 24
set output "vertex.eps"
set xlabel "{/Symbol e}_{/Times-Italic n}"
set ylabel "Re {/Symbol g}(i{/Symbol e}_{/Times-Italic n})"
set xrange[*:*]
set yrange[*:*]

set key left bottom

plot "<awk '$2 == 4 && $3 == 0' vertex_lo.dat"u 4:(($7-$9-$11+$13)/2.) title"{/Symbol w} = 0 (sp-zz)" w lp ls 1\
, "<awk '$2 == 4 && $3 == 0' vertex_tr.dat"u 4:9 title"(sp-pm)" w lp ls 2\
, "<awk '$2 == 4 && $3 == 5' vertex_lo.dat"u 4:(($7-$9-$11+$13)/2.) title"10{/Symbol p}{/Times-Italic T} (sp-zz)" w lp ls 3\
, "<awk '$2 == 4 && $3 == 5' vertex_tr.dat"u 4:9 title"(sp-pm)" w lp ls 4\

set output
