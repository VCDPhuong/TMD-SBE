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
plot "OUTPUT/Elec.txt" using 1:2 title 'valence' with line,\
     # "OUTPUT/Energy.txt" using 1:3 title 'conduct' with line
!set samples 40000
