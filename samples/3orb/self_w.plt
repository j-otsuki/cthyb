reset
set key bottom
set zeroaxis lw 2

# load a file in the script directory
SCRIPT_DIR = system("dirname ".ARG0)."/"
load SCRIPT_DIR."styles"


set terminal postscript eps color enhanced "Times-Roman" 24

set xlabel "{/Symbol-Oblique w}_{/Times-Italic n}"
set xrange[0:10]


stats "self_w.dat" u 1 nooutput
N = STATS_columns


set output "self_w_re.eps"
set ylabel "Re {/Symbol-Oblique S}_{imp}(i{/Symbol-Oblique w}_{/Times-Italic n})"
set yrange[*:*]

plot for [i=1:(N-1)/2] "self_w.dat"u 1:(column(2*i)) title sprintf("Re %d", i) w lp ls i



set output "self_w_im.eps"
set ylabel "Im {/Symbol-Oblique S}_{imp}(i{/Symbol-Oblique w}_{/Times-Italic n})"
set yrange[*:0]

plot for [i=1:(N-1)/2] "self_w.dat"u 1:(column(2*i+1)) title sprintf("Im %d", i) w lp ls i


set output
