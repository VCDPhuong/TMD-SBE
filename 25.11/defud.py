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
for nk1 in range(nkamax+1):
    for nk2 in range(nkamax+1):
        k=grid[nk1,nk2]
for nt in range(0,int((tmax-tmin)/deltat)):
    Elec[nt,0]=E0*exp(-(nt*deltat+tmin)**2/TauL**2)*cos(omegaL*(nt*deltat+tmin))*cos(phi)
    Elec[nt,1]=E0*exp(-(nt*deltat+tmin)**2/TauL**2)*cos(omegaL*(nt*deltat+tmin))*sin(phi)
    if nt !=0:
        Alec[nt,:]=Alec[nt-1,:]-Elec[nt,:]*deltat
def Elect(nt):
    return Elec[int(nt),:]+(Elec[int(nt)+1,:]-Elec[int(nt),:])*(nt-int(nt))
def Alect(nt):
    return Alec[int(nt),:]-Elec[int(nt),:]*(nt-int(nt))*deltat
def RHSsbeVG(R,nt):
    dummy=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
    colidee=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
    dummy=np.dot(p,Alect(nt))
    colidee=-np.einsum('ij...,ij->ij...',R,(1-np.identity(Norbital)),optimize=True)/T2
    return -1j/hbar*np.multiply(deltaE,R)-1j*qe/(hbar*me)*(np.einsum('ij...,jk...->ik...',dummy,R,optimize=True)-np.einsum('ij...,jk...->ik...',R,dummy,optimize=True))+colidee
def RK(rho,nt):
    RKdummy=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,4),dtype='complex128')
    RKdummy[...,0]=RHSsbeVG(rho,nt)
    RKdummy[...,1]=RHSsbeVG(rho+RKdummy[...,0]*deltat/2,nt+1/2)
    RKdummy[...,2]=RHSsbeVG(rho+RKdummy[...,1]*deltat/2,nt+1/2)
    RKdummy[...,3]=RHSsbeVG(rho+RKdummy[...,2]*deltat,nt+1)
    rho+=deltat/6*(RKdummy[...,0]+2*RKdummy[...,1]+2*RKdummy[...,2]+RKdummy[...,3])
    if np.any(np.diagonal(rho,axis1=0, axis2=1)>1.02):
        print("Density can't bigger than 1")
        exit()
    elif np.any(np.diagonal(rho,axis1=0, axis2=1)<-0.02):
        print("Density can't be negative")
        exit()
    return rho
def sbesolver():
    rho=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype="complex128")
    print("SBE. Solving...")
    for nt in range(0,int((tmax-tmin)/deltat)-1):
        start=time.time()
        if nt==0:
            with open("output/currentdensity.txt",mode ="a") as file:
                file.write(str(("-----------------------------------------------------------------------------------------------------"))+'\n')
            with open("output/density.txt",mode ="a") as file:
                file.write(str(("-----------------------------------------------------------------------------------------------------"))+'\n')
            if Norbital==3:
                rho[0,0,:,:]=1
            if Norbital==6:
                rho[0,0,:,:]=1
                rho[3,3,:,:]=1
        else:
            rho=RK(rho,nt)
        end=time.time()
        if np.divmod(nt,10/deltat)[1]==0: #report every 10 fs
            print("Time to calculate:", deltat,"[fs] is:","%.2f" % (end-start),"real time:","%.2f" % (nt*deltat+tmin),"on",tmax)
            print("Estimate to end:", time.ctime(end+(end-start)*(int((tmax-tmin)/deltat)-nt)))
            print(round(np.sum(rho[0,0,:,:]),2),round(np.sum(rho[1,1,:,:]),2),round(np.sum(rho[2,2,:,:]),2),round(np.sum(rho[3,3,:,:]),2),round(np.sum(rho[4,4,:,:]),2),round(np.sum(rho[5,5,:,:]),2))
            print("---------------------------------------------------------------------------------------------------------------------")
        if np.divmod(nt,0.2/deltat)[1]==0: #report current density every 1 fs
            Jdummykx=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
            Jdummyky=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
            J=np.zeros((2,2),dtype='complex128') #first index= 0=intra,1=inter; second:0 :kx, 1: ky
            Jdummykx=qe*np.multiply((p[...,0]-qe*Alec[nt,0]),np.transpose(rho,(1,0,2,3)))/me
            Jdummyky=qe*np.multiply((p[...,1]-qe*Alec[nt,1]),np.transpose(rho,(1,0,2,3)))/me
            J[0,0]=np.sum(np.trace(Jdummykx,axis1=0, axis2=1))
            J[0,1]=np.sum(np.trace(Jdummyky,axis1=0, axis2=1))
            J[1,0]=np.sum(Jdummykx)-np.sum(np.trace(Jdummykx,axis1=0, axis2=1))
            J[1,1]=np.sum(Jdummyky)-np.sum(np.trace(Jdummyky,axis1=0, axis2=1))
            with open("output/currentdensity.txt",mode ="a") as file:
                file.write(str(("-",np.round(nt*deltat+tmin,2),np.round(J[0,0],5),np.round(J[0,1],5),np.round(J[1,0],5),np.round(J[1,1],5),"-"))+'\n')
            with open("output/density.txt",mode ="a") as file:
                file.write(str(("-",np.round(nt*deltat+tmin,2),np.round(np.sum(rho[0,0,:,:]),5),np.round(np.sum(rho[1,1,:,:]),5),np.round(np.sum(rho[2,2,:,:]),5),np.round(np.sum(rho[3,3,:,:]),5),np.round(np.sum(rho[4,4,:,:]),5),np.round(np.sum(rho[5,5,:,:]),5),"-"))+'\n')
    print("SBE Calculation: Done")