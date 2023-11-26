import sys
sys.path.append('Z:/root/TB')
from Hamiltonian import H,Vpkx,Vpky,Norbital,a
from variable import variable,nkamax
import Hamiltonian
import numpy as np
from numpy import linalg as LA
from numpy import sin, cos, sqrt, pi, exp
import matplotlib.pyplot as plt
import time
from numba import jit,cuda
import multiprocessing as mp
from multiprocessing import Process, Pipe, Pool
import csv
#time in fs
tmin=-300
tmax=1000
deltat=0.02 #fs
tpeak=0
#Elect Interaction
E0=1e-2 #V/nm
phi=0 #rad
TauL=60 #fs
hbar=0.6582 #[eV.fs]
epsilonL=1.8/4 #[eV]
omegaL=epsilonL/hbar # [1/fs]
qe=-1 #[unit charge]
me=9.1094/1.6022 #[eV.fs^2/nm^2]; E=mc^2;[kg]=[J(s/m)^2]=1/qe [eV (s/m)^2]
#fourier transform for harmonic
NomegaL=40
deltaomega=0.05
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
for nk1 in range(nkamax+1):
    for nk2 in range(nkamax+1):
        k=grid[nk1,nk2]
        E[:,nk1,nk2],V[:,:,nk1,nk2]=LA.eigh(H(k[0],k[1]))
for nt in range(0,int((tmax-tmin)/deltat)):
    Elec[nt,0]=E0*exp(-(nt*deltat+tmin-tpeak)**2/TauL**2)*cos(omegaL*(nt*deltat+tmin-tpeak))*cos(phi)
    Elec[nt,1]=E0*exp(-(nt*deltat+tmin-tpeak)**2/TauL**2)*cos(omegaL*(nt*deltat+tmin-tpeak))*sin(phi)
    if nt !=0:
        Alec[nt,:]=Alec[nt-1,:]-Elec[nt,:]*deltat
def Alect(nt):
    return Alec[int(nt),:]-Elec[int(nt),:]*(nt-int(nt))*deltat
def RHSsbeVG(R,nt):
    dummy=p[...,:]@Alect(nt)
    RHS=-1j/hbar*np.multiply(deltaE,R)-1j*qe/(hbar*me)*(np.einsum("ij...,jk...->ik...",dummy,R)-np.einsum("ij...,jk...->ik...",R,dummy))
    return RHS
def RK(rho,nt):
    RKdummy=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,4),dtype='complex128')
    RKdummy[...,0]=RHSsbeVG(rho,nt)
    RKdummy[...,1]=RHSsbeVG(rho+RKdummy[...,0]*deltat/3,nt+1/3)
    RKdummy[...,2]=RHSsbeVG(rho-RKdummy[...,0]*deltat/3+RKdummy[...,1]*deltat,nt+2/3)
    RKdummy[...,3]=RHSsbeVG(rho+RKdummy[...,0]*deltat-RKdummy[...,1]*deltat+RKdummy[...,2]*deltat,nt+1)
    rho+=deltat/8*(RKdummy[...,0]+3*RKdummy[...,1]+3*RKdummy[...,2]+RKdummy[...,3])
    rhodiagcheck=np.diagonal(rho, axis1=0, axis2=1)
    if np.any(np.real(rhodiagcheck)>1.02):
        print("Density can't bigger than 1")
        exit()
    elif np.any(np.real(rhodiagcheck)<-2e-2):
        print("Density can't be negative")
        exit()
    return rho
def sbesolver():
    rho=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype="complex128")
    Jdummy=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex128')
    Jdummyintra=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex128')
    Jdummyinter=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex128')
    print("SBE. Solving...")
    for nt in range(0,int((tmax-tmin)/deltat)):
        start=time.time()
        if nt==0:
            if Norbital==3:
                rho[0,0,:,:]=1
            if Norbital==6:
                rho[0,0,:,:]=1
                rho[1,1,:,:]=1
        else:
            rho=RK(rho,nt)
        end=time.time()
        if np.divmod(nt,10/deltat)[1]==0: #report every 10 fs
            print("Time to calculate:", deltat,"[fs] is:","%.2f" % (end-start),"real time:","%.2f" % (nt*deltat+tmin),"on",tmax)
            print("Estimate to end:", time.ctime(end+(end-start)*(int((tmax-tmin)/deltat)-nt)))
            print(round(np.sum(rho[0,0,:,:]),2),round(np.sum(rho[1,1,:,:]),2),round(np.sum(rho[2,2,:,:]),2),round(np.sum(rho[3,3,:,:]),2),round(np.sum(rho[4,4,:,:]),2),round(np.sum(rho[5,5,:,:]),2))
            print("-------------------------------------------------------------------------------------------------")
            Jdummy=qe/me*np.einsum('ijkl...,ijkl->ijkl...',p[...,:]-qe*Alec[nt,:],rho)
            Jdummyintra=np.einsum('ij...,ij->ij...',Jdummy,Delta) #(einsum inside give (Nor,Nor,nkamax+1,nkamax+1,2))
            Jdummyinter=np.einsum('ij...,ij->ij...',Jdummy,1-Delta) #(einsum inside give (Nor,Nor,nkamax+1,nkamax+1,2))
            J=np.zeros((2,2),dtype='complex128') #first index= 0=intra,1=inter; second:0 :kx, 1: ky
            J[0,0]=np.sum(Jdummyintra[...,0])
            J[0,1]=np.sum(Jdummyintra[...,1])
            J[1,0]=np.sum(Jdummyinter[...,0])
            J[1,1]=np.sum(Jdummyinter[...,1])
            with open("output/currentdensity.csv",mode ="+a") as write:
                writercur = csv.writer(open('output/Currentdensity.csv', 'a', newline =''), quoting = csv.QUOTE_ALL,delimiter =';')
                writercur.writerows(np.transpose([[nt*deltat+tmin],[J[0,0]],[J[0,1]],[J[1,0]],[J[1,1]]]))