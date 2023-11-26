# Gnuplot script file for plotting data in file "force.dat"
# This file is called   force.p
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Force Deflection Data for a Beam and a Column"
set xlabel "Step calculate"
set ylabel "value"
set grid
plot "Energyonkx.txt" using 1:3 title 'VB1' with points,\
    "Energyonkx.txt" using 1:4 title 'VB2' with points,\
    "Energyonkx.txt" using 1:5 title 'CB1' with points,\
    "Energyonkx.txt" using 1:6 title 'CB2' with points,\
    "Energyonkx.txt" using 1:7 title 'CB3' with points,\
    "Energyonkx.txt" using 1:8 title 'CB4' with points
set samples 40000