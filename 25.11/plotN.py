import numpy as np
import matplotlib.pyplot as plt
#in file col 1 is time, col 2 is intra kx, 3 is intra ky, 4 is inter kx, 5 is inter ky, 0 and 6 is dummy
t=  np.loadtxt("output/density.txt",usecols=1,skiprows=1,delimiter=",",dtype='complex')
Nv1=np.loadtxt("output/density.txt",usecols=2,skiprows=1,delimiter=",",dtype='complex')
Nc1=np.loadtxt("output/density.txt",usecols=3,skiprows=1,delimiter=",",dtype='complex')
Nc2=np.loadtxt("output/density.txt",usecols=4,skiprows=1,delimiter=",",dtype='complex')
Nv2=np.loadtxt("output/density.txt",usecols=5,skiprows=1,delimiter=",",dtype='complex')
Nc3=np.loadtxt("output/density.txt",usecols=6,skiprows=1,delimiter=",",dtype='complex')
Nc4=np.loadtxt("output/density.txt",usecols=7,skiprows=1,delimiter=",",dtype='complex')
Nc=Nc1+Nc2+Nc3+Nc4
Nv=Nv1+Nv2
plt.plot(t,Nc,'r')
#plt.plot(t,Nv,'xb')
#plt.plot(Jffreq,I[:,2],'s')
plt.grid()
plt.show()