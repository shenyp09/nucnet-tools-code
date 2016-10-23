set term pdf font 'Times-Roman'
set output "al26.pdf"

unset key
set format y "10^{%T}"
set ylabel "" font 'Times-Roman,12'   offset -3.0,0
# set xlabel "log_{10}(factor)" font 'Times-Roman, 12' offset 0,-0.5
set xtics font 'Times-Roman,9'
set ytics font 'Times-Roman,9'
plot 'data.txt' u (log10($1)):2 smooth csplines lc 1, '' u (log10($1)):2 w p pt 7 pointsize 1 lc 1
