reset
set key bottom
set zeroaxis lw 2

# load a file in the script directory
SCRIPT_DIR = system("dirname ".ARG0)."/"
load SCRIPT_DIR."styles"


set terminal postscript eps color enhanced "Times-Roman" 24
set output "stat.eps"
set xlabel "Expansion order"
set ylabel ""
set xrange[*:*]
set yrange[0:*]

stats "stat.dat" u 1 nooutput
N = STATS_columns

plot for [i=1:N-1] "stat.dat"u 1:(column(i+1)) title sprintf("%d", i) w lp ls i

set output
