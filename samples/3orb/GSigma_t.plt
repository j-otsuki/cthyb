reset
set key bottom
set zeroaxis lw 2

# load a file in the script directory
SCRIPT_DIR = system("dirname ".ARG0)."/"
load SCRIPT_DIR."styles"


set terminal postscript eps color enhanced "Times-Roman" 24
set output "GSigma_t.eps"
set xlabel "{/Symbol-Oblique t}"
set ylabel "({/Times-Italic G}{/Symbol-Oblique S})({/Symbol-Oblique t})"
set xrange[0:*]
set yrange[:0]

stats "GSigma_t.dat" u 1 nooutput
N = STATS_columns

plot for [i=1:(N-1)/2] "GSigma_t.dat"u 1:(column(2*i)):(column(2*i+1)) title sprintf("%d", i) w ye ls i

set output
