import numpy as np
import matplotlib.pyplot as plt
from defud import hbar,qe,me,tmin,tmax
epsilonL=1.81
omegaL=epsilonL/(hbar)
NomegaL=30
deltaomega=0.1
deltat=0.2
from numpy import sin, cos, sqrt, pi, exp
#time domaino
tdomainmin=-200
ntdmin=int(tdomainmin/deltat)
tdomainmax=200
ntdmax=int((tdomainmax+tdomainmin)/deltat)
#in file col 1 is time, col 2 is intra kx, 3 is intra ky, 4 is inter kx, 5 is inter ky, 0 and 6 is dummy
t=       np.loadtxt("output/currentdensity.txt",usecols=1,skiprows=1,delimiter=",",dtype='complex')
Jintrakx=np.loadtxt("output/currentdensity.txt",usecols=2,skiprows=1,delimiter=",",dtype='complex')
Jintraky=np.loadtxt("output/currentdensity.txt",usecols=3,skiprows=1,delimiter=",",dtype='complex')
Jinterkx=np.loadtxt("output/currentdensity.txt",usecols=4,skiprows=1,delimiter=",",dtype='complex')
Jinterky=np.loadtxt("output/currentdensity.txt",usecols=5,skiprows=1,delimiter=",",dtype='complex')
for nt in range(len(t)):
    if t[nt]<tdomainmin or t[nt]>tdomainmax:
        Jintrakx[nt]=0
        Jintraky[nt]=0
        Jinterkx[nt]=0
        Jinterky[nt]=0
Jtotalkx=Jintrakx+Jinterkx
Jtotalky=Jintraky+Jinterky
I=np.zeros((int(NomegaL/deltaomega),6),dtype='complex128')
Jffreq=np.zeros((int(NomegaL/deltaomega)),dtype='float32')
print("Elec calculation: Done")
for i in range(0,int(NomegaL/deltaomega)):
    Jf=np.zeros((6),dtype='complex128')
    for nt in range(0,len(Jtotalkx)):
        Jf[0]+=Jintrakx[nt]*exp(-1j*i*deltaomega*omegaL*(t[nt]))*deltat
        Jf[1]+=Jintraky[nt]*exp(-1j*i*deltaomega*omegaL*(t[nt]))*deltat
        Jf[2]+=Jinterkx[nt]*exp(-1j*i*deltaomega*omegaL*(t[nt]))*deltat
        Jf[3]+=Jinterky[nt]*exp(-1j*i*deltaomega*omegaL*(t[nt]))*deltat
        Jf[4]+=Jtotalkx[nt]*exp(-1j*i*deltaomega*omegaL*(t[nt]))*deltat
        Jf[5]+=Jtotalky[nt]*exp(-1j*i*deltaomega*omegaL*(t[nt]))*deltat
    I[i,:]=np.abs(Jf[:])**2
    Jffreq[i]=i*deltaomega
    print(i,"on", int(NomegaL/deltaomega))
plt.title("Fourier of CurrentDensity")
plt.yscale("log")
plt.xlabel(r'frequency [$\omega_L$]')
plt.xlabel(r'Intensity [$a.u$]')
plt.plot(Jffreq,I[:,0],'b',label='Intra Current')
plt.plot(Jffreq,I[:,2],'r',label='Inter Current')
plt.plot(Jffreq,I[:,4],'y',label='Total Current')
plt.legend()
plt.grid()
#plt.savefig('currentdensity(0.326,1e0).png')
plt.show()