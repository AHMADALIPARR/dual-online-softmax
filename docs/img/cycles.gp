# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
# Copyright Ahmad Ali Parr
# Cycles per tile vs keys per row, 4 rows (BR=4), BC=4 keys per tile.
# Data: cycles.dat, written by sim/fa_dual/tb_fa_dual.sv (make figs copies it here). Render: gnuplot cycles.gp, or gnuplot -e "fmt='png'" cycles.gp
if (!exists("fmt")) fmt = "svg"
if (fmt eq "png") { set terminal pngcairo size 900,520 font "Helvetica,12"; set output "cycles.png" } else { set terminal svg size 900,520 font "Helvetica,12"; set output "cycles.svg" }
set title "fa\\_online\\_softmax vs fa\\_dual\\_online\\_softmax (BR=4 rows, BC=4 keys/tile)"
set xlabel "keys per row"
set ylabel "cycles per tile = total / (rows * keys / 4)"
set logscale x 2
set yrange [0:13]
set key top right
set grid
f(c,r,k) = c / (r * k / 4.0)
mod = "< awk '$1==4 && $3==0' cycles.dat"
ris = "< awk '$1==4 && $3==1' cycles.dat"
plot mod using 2:(f($4,$1,$2)) with linespoints lw 2 title "original, contract issue", \
     mod using 2:(f($5,$1,$2)) with linespoints lw 2 title "original, safe pipelined issue (BC+1)", \
     mod using 2:(f($6,$1,$2)) with linespoints lw 2 title "dual, ping-pong, moderate logits", \
     ris using 2:(f($6,$1,$2)) with linespoints lw 2 dt 2 title "dual, ping-pong, rising logits (max moves every tile)"
