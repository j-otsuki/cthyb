reset
set key
set zeroaxis lw 3

set style line 1 lt 1 lw 2 pt 1 ps 1.5
set style line 2 lt 2 lw 2 pt 2 ps 1.5
set style line 3 lt 3 lw 2 pt 6 ps 1.5
set style line 4 lt 4 lw 2 pt 4 ps 1.5


set terminal postscript eps color enhanced "Times-Roman" 24
set output "chi_t.eps"
set xlabel "{/Symbol t}"
set ylabel "{/Symbol c}_{imp}({/Symbol t}) / {/Times-Italic C_N}"
set xrange[0:]
set yrange[:]

plot "chi_t.dat"u 1:2:3 title"sp" w ye ls 1\
, ""u 1:4:5 title"ch" w ye ls 2\
, ""u 1:6:7 title"00" w ye ls 3\
, ""u 1:8:9 title"01" w ye ls 4\

set output


set terminal postscript eps color enhanced "Times-Roman" 24
set output "chi_w.eps"
set xlabel "{/Symbol w}_{/Times-Italic n}"
set ylabel "Re {/Symbol c}_{imp}(i{/Symbol w}_{/Times-Italic n}) / {/Times-Italic C_N}"
set xrange[0:1]
set yrange[*:*]

plot "chi_w.dat"u 1:2 title"sp" w lp ls 1\
, "chi_w.dat"u 1:3 title"ch" w lp ls 2\
, "chi_w.dat"u 1:4 title"00 (re)" w lp ls 3\
, "chi_w.dat"u 1:6 title"01 (re)" w lp ls 4\

set output


# set terminal postscript eps color enhanced "Times-Roman" 24
# set output "chi_pade.eps"
# set xlabel "{/Symbol w}"
# set ylabel "{/Symbol c}_{imp}({/Symbol w}) / {/Times-Italic C_N}"
# set xrange[0:0.5]
# set yrange[*:*]

# plot "chi_pade.dat"u 1:3 title"sp" w l ls 1\
# , ""u 1:5 title"ch" w l ls 2\

# set output


# set terminal postscript eps color enhanced "Times-Roman" 24
# set output "chi_pade2.eps"
# set xlabel "{/Symbol w}"
# set ylabel "Im {/Symbol c}_{imp}({/Symbol w}) / {/Symbol w} {/Times-Italic C_N}"
# set yrange[*:*]

# plot "chi_pade.dat"u 1:($3/$1) title"sp" w l ls 1\
# , ""u 1:($5/$1) title"ch" w l ls 2\

# set output
