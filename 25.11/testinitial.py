import sys
#sys.path.append('Z:/root/TB')
from Hamiltonian import H,Vpkx,Vpky,Norbital,a
from variable import variable,nkamax
import numpy as np
from numpy import linalg as LA
from numpy import sin, cos, sqrt, pi, exp
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})
import time
import csv
#time in fs
tmin=-400
tmax=400
deltat=0.02 #fs
tpeak=0
#Elect Interaction
E0=1e0 #V/nm
phi=0 #rad
TauL=60 #fs
hbar=0.6582 #[eV.fs]
epsilonL=1.8/20 #[eV]
omegaL=epsilonL/(hbar)# [1/fs]
qe=-1 #[unit charge]
me=9.1094/1.6022 #[eV.fs^2/nm^2]; E=mc^2;[kg]=[J(s/m)^2]=1/qe [eV (s/m)^2]
#fourier transform for harmonic
NomegaL=40
deltaomega=0.05
Elec=np.zeros((int((tmax-tmin)/deltat),2),dtype='float32')
Alec=np.zeros((int((tmax-tmin)/deltat),2),dtype='float32')
for nt in range(0,int((tmax-tmin)/deltat)):
    Elec[nt,0]=E0*exp(-(nt*deltat+tmin-tpeak)**2/TauL**2)*cos(omegaL*(nt*deltat+tmin-tpeak))*cos(phi)
    Elec[nt,1]=E0*exp(-(nt*deltat+tmin-tpeak)**2/TauL**2)*cos(omegaL*(nt*deltat+tmin-tpeak))*sin(phi)
    if nt !=0:
        Alec[nt,:]=Alec[nt-1,:]-Elec[nt,:]*deltat
def Alect(nt):
    return Alec[int(nt),:]-Elec[int(nt),:]*(nt-int(nt))*deltat
p=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex64')
Xi=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex64')
V=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex64')
rho=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
RKdummy=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,4),dtype='complex128')
N=np.zeros((Norbital,int((tmax-tmin)/deltat)),dtype='float32')
J=np.zeros((int((tmax-tmin)/deltat),2,3),dtype='complex128') #second index: 0 is x,1 is y. Final Index:0 is intra, 1 is inter, 3 is total
BerryCurvz=np.zeros((Norbital,nkamax+1,nkamax+1),dtype='float64')
Delta=np.zeros((Norbital,Norbital),dtype='int16')
deltaE=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='float32')
E=np.zeros((Norbital,nkamax+1,nkamax+1),dtype='float64')
grid,p,Xi,BerryCurvz,Delta,deltaE,E,V=variable()
dummy1=np.zeros((Norbital,nkamax+1),dtype='complex128')
dummy2=np.zeros((nkamax+1,2),dtype='complex128')
#RHS=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
for nk in range(nkamax+1):
    dummy2[nk,:]=grid[nk,nk,:]
    for nu in range(Norbital):
        dummy1[nu,nk]=BerryCurvz[nu,nk,nk]
plt.plot(dummy2[:,0],dummy1[0,:],'xb')
#plt.plot(grid[:,:,0],np.transpose(grid[:,:,1],(1,0)),'b')
plt.xticks([dummy2[0,0],dummy2[int(nkamax/6)+1,0],dummy2[int(nkamax/2)+1,0],dummy2[int(5*nkamax/6),0],dummy2[nkamax,0]],["M","-K",r"$\Gamma$","K","M"])
plt.grid()
plt.ylabel("Energy [eV]")
plt.show()