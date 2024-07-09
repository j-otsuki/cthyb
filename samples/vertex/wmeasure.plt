reset
set key
set zeroaxis lw 2

set style line 1 lt 1 lw 2 pt 1 ps 1.5
set style line 2 lt 2 lw 2 pt 2 ps 1.5
set style line 3 lt 3 lw 2 pt 6 ps 1.5
set style line 4 lt 4 lw 2 pt 4 ps 1.5



set terminal postscript eps color enhanced "Times-Roman" 24
set output "wmeasure_Gf_w.eps"
set xlabel "{/Symbol e}_{/Times-Italic n}"
set ylabel "{/Times-Italic G}_{imp}(i{/Symbol e}_{/Times-Italic n})"
set xrange[0:10]
set yrange[*:*]

set key bottom

plot "Gf_w.dat"u 1:2 title"0 Re" ls 1, ""u 1:3 title"0 Im" ls 2\
, "wmeasure_Gf_w.dat"u 1:2 title"(w-measure) 0 Re" ls 3, ""u 1:3 title"0 Im" ls 4\

# , ""u 1:4 title"1 Re" ls 3, ""u 1:5 title"1 Im" ls 4\
# , ""u 1:4 title"1 Re" ls 3, ""u 1:5 title"1 Im" ls 4\

set output




set terminal postscript eps color enhanced "Times-Roman" 24
set output "wmeasure_self_w.eps"
set xlabel "{/Symbol e}_{/Times-Italic n}"
set ylabel "{/Symbol S}_{imp}(i{/Symbol e}_{/Times-Italic n})"
set xrange[0:10]
set yrange[*:*]

set key bottom

plot "self_w.dat"u 1:2 title"0 Re" ls 1, ""u 1:3 title"0 Im" ls 2\
, "wmeasure_self_w.dat"u 1:2 title"(w-measure) 0 Re" ls 3, ""u 1:3 title"0 Im" ls 4\

set output





set terminal postscript eps color enhanced "Times-Roman" 24
set output "wmeasure_chi_w.eps"
set xlabel "{/Symbol w}_{/Times-Italic n}"
set ylabel "{/Symbol c}_{imp}(i{/Symbol w}_{/Times-Italic n})"
set xrange[0:10]
set yrange[*:*]

set key bottom

plot "chi_w.dat"u 1:4 title"00 Re" ls 1, ""u 1:6 title"01 Re" ls 2\
, "wmeasure_chi_w.dat"u 1:2 title"(w-measure) 00 Re" ls 3, ""u 1:4 title"01 Re" ls 4\

set output
