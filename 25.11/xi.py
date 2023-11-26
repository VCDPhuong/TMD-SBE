import sys
sys.path.append('Z:/root/TB')
from variableud import variable,nkamax
from hamud import Hamu,Hamd,Vpkx,Vpky,Norbital
import numpy as np
from numpy import linalg as LA
from numpy import sin, cos, sqrt, pi, exp
import matplotlib.pyplot as plt
import time
import csv
#time in fs
tmin=-400
tmax=400
deltat=0.02 #fs
#Elect Interaction
E0=1e-2 #V/nm
phi=0 #rad
theta=0 #rad
TauL=60*sqrt(2) #fs (the TauL is not have the form x^2/2\tauL**2)
hbar=0.6582 #[eV.fs]
epsilonL=1.81 #[eV]
omegaL=epsilonL/(hbar)# [1/fs]
qe=1 #[unit charge]
me=9.1094/1.6022 #[eV.fs^2/nm^2]; E=mc^2;[kg]=[J(s/m)^2]=1/qe [eV (s/m)^2]
#dephasing
T2=5.2 #[fs]
#fourier transform for harmonic
Elec=np.zeros((int((tmax-tmin)/deltat),2),dtype='float32')
Alec=np.zeros((int((tmax-tmin)/deltat),2),dtype='float32')
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
print("Variable: Defined")
if Norbital==3:
    valence=[0]
    conduct=[1,2]
elif Norbital==6:
    valence=[0,3]
    conduct=[1,2,4,6]
print("Valance,conduct index: defined")
for nk1 in range(nkamax+1):
    for nk2 in range(nkamax+1):
        k=grid[nk1,nk2]
for nt in range(0,int((tmax-tmin)/deltat)):
    Elec[nt,0]=E0*exp(-(nt*deltat+tmin)**2/TauL**2)*cos(omegaL*(nt*deltat+tmin))*cos(phi)*cos(theta)
    Elec[nt,1]=E0*exp(-(nt*deltat+tmin)**2/TauL**2)*cos(omegaL*(nt*deltat+tmin))*sin(phi)*cos(theta)
    if nt !=0:
        Alec[nt,:]=Alec[nt-1,:]-Elec[nt,:]*deltat
def Elect(nt):
    return Elec[int(nt),:]+(Elec[int(nt)+1,:]-Elec[int(nt),:])*(nt-int(nt))
def Alect(nt):
    return Alec[int(nt),:]-Elec[int(nt),:]*(nt-int(nt))*deltat
