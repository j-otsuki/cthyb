reset
set key top
set zeroaxis lw 2

# load a file in the script directory
SCRIPT_DIR = system("dirname ".ARG0)."/"
load SCRIPT_DIR."styles"


set terminal postscript eps color enhanced "Times-Roman" 24
set output "Gf_pade.eps"
set xlabel "{/Symbol-Oblique w}"
set ylabel "{/Symbol-Oblique r}_{imp}({/Symbol-Oblique w})"
set xrange[*:*]
set yrange[0:*]

stats "Gf_pade.dat" u 1 nooutput
N = STATS_columns

plot for [i=1:(N-1)/2] "Gf_pade.dat"u 1:(-column(2*i+1)/pi) title sprintf("%d", i) w l ls i

set output
