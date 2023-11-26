import numpy as np
import matplotlib.pyplot as plt
from defud import hbar,qe,me,tmin,tmax,Alect, deltat,Elect
from numpy import sin, cos, sqrt, pi, exp
#time domaino
#in file col 1 is time, col 2 is intra kx, 3 is intra ky, 4 is inter kx, 5 is inter ky, 0 and 6 is dummy
t=       np.loadtxt("output/currentdensity(1.81,1e-2).txt",usecols=1,skiprows=1,delimiter=",",dtype='complex')
Jintrakx=np.loadtxt("output/currentdensity(1.81,1e-2).txt",usecols=2,skiprows=1,delimiter=",",dtype='complex')
Jintraky=np.loadtxt("output/currentdensity(1.81,1e-2).txt",usecols=3,skiprows=1,delimiter=",",dtype='complex')
Jinterkx=np.loadtxt("output/currentdensity(1.81,1e-2).txt",usecols=4,skiprows=1,delimiter=",",dtype='complex')
Jinterky=np.loadtxt("output/currentdensity(1.81,1e-2).txt",usecols=5,skiprows=1,delimiter=",",dtype='complex')
Jtotalkx=Jintrakx+Jinterkx
Jtotalky=Jintraky+Jinterky
fig, (ax1, ax2) = plt.subplots(2, sharex=True)
for nt in range(len(t)):
    ax1.plot(t[nt], Alect((t[nt]-tmin)/deltat)[0],'bo')
ax1.set_title("Potential of Electrodynamic wave")
ax1.grid()
ax2.plot(t, Jintrakx,'b--',label='Intra Current')
ax2.plot(t, Jinterkx,'r--',label='Inter Current')
ax2.plot(t, Jtotalkx,'y--',label='Total Current')
ax2.legend()
ax2.grid()
ax2.set_title("Current Density")
plt.show()