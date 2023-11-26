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
plot "px.txt" using 1:3 title 'px00' with points,\
     "px.txt"    using 1:4 title 'px01' with points,\
     "px.txt"    using 1:5 title 'px02' with points
#set samples 40000