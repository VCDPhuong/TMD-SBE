import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.append('Z:/root/TB')
from defud import tmin,tmax,E0,phi,TauL,omegaL,qe,me
NomegaL=30
deltaomega=0.1
deltat=0.2 #txt file is report every 0.2 fs
from numpy import exp
#in file col 1 is time, col 2 is intra kx, 3 is intra ky, 4 is inter kx, 5 is inter ky, 0 and 6 is dummy
t=np.loadtxt("output/currentdensity(1.81,1e0).txt",usecols=1,skiprows=1,delimiter=",",dtype='complex')
Jintrakx=np.loadtxt("output/currentdensity(1.81,1e0).txt",usecols=2,skiprows=1,delimiter=",",dtype='complex')
Jintraky=np.loadtxt("output/currentdensity(1.81,1e0).txt",usecols=3,skiprows=1,delimiter=",",dtype='complex')
Jinterkx=np.loadtxt("output/currentdensity(1.81,1e0).txt",usecols=4,skiprows=1,delimiter=",",dtype='complex')
Jinterky=np.loadtxt("output/currentdensity(1.81,1e0).txt",usecols=5,skiprows=1,delimiter=",",dtype='complex')
Jtotalkx=Jintrakx+Jinterkx
Jtotalky=Jintraky+Jinterky
Jf=np.zeros((int(NomegaL/deltaomega),4),dtype='complex128')
I=np.zeros((int(NomegaL/deltaomega),4),dtype='complex128')
Jffreq=np.linspace(0,NomegaL,num=int(NomegaL/deltaomega))
print("Elec calculation: Done")
for i in range(0,int(NomegaL/deltaomega)):
    for nt in range(0,len(Jtotalkx)):
        Jf[i,0]+=Jintrakx[nt]*exp(-1j*i*deltaomega*omegaL*(t[nu]))*deltat
        Jf[i,1]+=Jintraky[nt]*exp(-1j*i*deltaomega*omegaL*(t[nu]))*deltat
        Jf[i,2]+=Jinterkx[nt]*exp(-1j*i*deltaomega*omegaL*(t[nu]))*deltat
        Jf[i,3]+=Jinterky[nt]*exp(-1j*i*deltaomega*omegaL*(t[nu]))*deltat
    I[i,:]=np.abs(Jf[i,:])**2
    Jffreq[i]=i*deltaomega
    print(i,"on", int(NomegaL/deltaomega))
plt.yscale("log")
plt.plot(Jffreq,I[:,0],'b')
plt.plot(Jffreq,I[:,2],'r')
plt.grid()
plt.show()