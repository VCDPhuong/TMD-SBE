set xrange[-200:1000]
set xlabel "time"
set multiplot layout 3,1 columns title "Bigg E"
set title 'Electrical Field'
plot "OUTPUT/Elec.txt" using 1:2 title 'Electrical Field' with line
set title 'Polarization'
plot "OUTPUT/Polarize.txt" using 1:2 title 'Polarization' with line
set title 'Density'
plot "OUTPUT/density.txt" using 1:4 title 'Density' with line
set grid