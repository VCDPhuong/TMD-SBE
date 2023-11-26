import numpy as np
import math
from scipy.linalg import svd
from numpy import linalg as LA
from numpy import sin, cos, sqrt, pi, exp
import matplotlib.pyplot as plt
import csv
import time
from numba import jit,cuda
a=31.9 #nanomet
# cac gia tri hoping
e1=0.683
e2=1.707
t0=-0.146
t1=-0.114
t2=0.506
t11=0.085
t12=0.162
t22=0.073
r0=0.060
r1=-0.236
r2=0.067
r11=0.016
r12=0.087
u0=-0.038
u1=0.046
u2=0.001
u11=0.266
u12=-0.176
u22=-0.150
lam=0.073
nkamax=49
#time in fs
tmin=-200
tmax=200
deltat=0.02 #fs
tpeak=0
#Elect Interaction
E0=1e0 #V/nm
phi=0 #rad
TauL=60 #fs
hbar=0.6582 #[eV.fs]
epsilonL=0.589 #[eV]
omegaL=epsilonL/(hbar)# [1/fs]
qe=-1 #[unit charge]
me=9.1094/1.6022 #[eV.fs^2/nm^2]; E=mc^2;[kg]=[J(s/m)^2]=1/qe [eV (s/m)^2]
#fourier transform for harmonic
NomegaL=40
deltaomega=0.05
#define electronmagnectic waves:
Elec=np.zeros((int((tmax-tmin)/deltat),2),dtype='float32')
Alec=np.zeros((int((tmax-tmin)/deltat),2),dtype='float32')
for nt in range(0,int((tmax-tmin)/deltat)):
    Elec[nt,0]=E0*exp(-(nt*deltat+tmin-tpeak)**2/TauL**2)*cos(omegaL*(nt*deltat+tmin-tpeak))*cos(phi)
    Elec[nt,1]=E0*exp(-(nt*deltat+tmin-tpeak)**2/TauL**2)*cos(omegaL*(nt*deltat+tmin-tpeak))*sin(phi)
    if nt !=0:
        Alec[nt,:]=Alec[nt-1,:]-Elec[nt,:]*deltat
#in case nt not a integer:
@jit(target_backend='cuda',nopython=True)
def Elect(t):
    return E0*exp(-(t+tmin-tpeak)**2/TauL**2)*cos(omegaL*(t+tmin-tpeak))*[[cos(phi)],[sin(phi)]]
@jit(target_backend='cuda',nopython=True)
def Alect(t):
    return Alec[nt-1,:]-Elec[nt,:]*(t-math.floor((t-tmin)/deltat)*deltat)
print("Electronmagnectic waves: defined")
#Hamiltonian:
@jit(nopython=True)
def H(kx,ky):
    h0=e1+2*t0*(2*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*kx*a/2))+2*r0*(2*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*ky*a*sqrt(3)/2))+2*u0*(2*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+cos(4*kx*a/2))
    h1=complex(-2*sqrt(3.)*t2*sin(kx*a/2)*sin(ky*a*sqrt(3)/2)+2*(r1+r2)*sin(3*kx*a/2)*sin(ky*a*sqrt(3)/2)-2*sqrt(3.)*u2*sin(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2) , 2*t1*sin(kx*a/2)*(2*cos(kx*a/2)+cos(ky*a*sqrt(3)/2))+2*(r1-r2)*sin(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+2*u1*sin(2*kx*a/2)*(2*cos(2*kx*a/2)+cos(2*ky*a*sqrt(3)/2)))
    h2=complex(2*t2*(cos(2*kx*a/2)-cos(kx*a/2)*cos(ky*a*sqrt(3)/2))-2/sqrt(3.)*(r1 + r2)*(cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)-cos(2*ky*a*sqrt(3)/2))+2*u2*(cos(4*kx*a/2)-cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)),2*sqrt(3.)*t1*cos(kx*a/2)*sin(ky*a*sqrt(3)/2)+2/sqrt(3.)*sin(ky*a*sqrt(3)/2)*(r1-r2)*(cos(3*kx*a/2)+2*cos(ky*a*sqrt(3)/2))+2*sqrt(3.)*u1*cos(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2))
    h11=e2+(t11+3*t22)*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+2*t11*cos(2*kx*a/2)+4*r11*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+2*(r11+sqrt(3.)*r12)*cos(2*ky*a*sqrt(3)/2)+(u11+3*u22)*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+2*u11*cos(4*kx*a/2)
    h12=complex(sqrt(3.)*(t22-t11)*sin(kx*a/2.)*sin(ky*a*sqrt(3)/2.)+4*r12*sin(3*kx*a/2)*sin(ky*a*sqrt(3)/2)+sqrt(3.)*(u22-u11)*sin(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2),4*t12*sin(kx*a/2)*(cos(kx*a/2)-cos(ky*a*sqrt(3)/2))+4*u12*sin(2*kx*a/2)*(cos(2*kx*a/2)-cos(2*ky*a*sqrt(3)/2)))
    h22=e2+(3*t11+t22)*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+2*t22*cos(2*kx*a/2)+2*r11*(2*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*ky*a*sqrt(3)/2))+2/sqrt(3.)*r12*(4*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)-cos(2*ky*a*sqrt(3)/2))+(3*u11+u22)*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+2*u22*cos(4*kx*a/2)
    return [[h0,h1,h2,0,0,0],[np.conj(h1),h11,h12+1j*lam,0,0,0],[np.conj(h2),np.conj(h12)-1j*lam,h22,0,0,0],[0,0,0,h0,h1,h2],[0,0,0,np.conj(h1),h11,h12-1j*lam],[0,0,0,np.conj(h2),np.conj(h12)+1j*lam,h22]]
Norbital,col=np.shape(H(0,0))
@jit(nopython=True)
def Vpkx(kx,ky):
    delta=1e-14
    return (np.array(H(kx+delta,ky), dtype='complex128')-np.array(H(kx-delta,ky),dtype='complex128'))/(2*delta)
@jit(nopython=True)
def Vpky(kx,ky):
    delta=1e-14
    return (np.array(H(kx,ky+delta), dtype='complex128')-np.array(H(kx,ky-delta),dtype='complex128'))/(2*delta)
p=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex64')
Xi=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex64')
V=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex64')
rho=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
RKdummy=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,4),dtype='complex128')
RHS=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
N=np.zeros((Norbital,int((tmax-tmin)/deltat)),dtype='float32')
J=np.zeros((int((tmax-tmin)/deltat),2,3),dtype='complex128') #second index: 0 is x,1 is y. Final Index:0 is intra, 1 is inter, 3 is total
BerryCurvz=np.zeros((Norbital,nkamax+1,nkamax+1),dtype='float64')
Delta=np.zeros((Norbital,Norbital),dtype='int16')
E=np.zeros((Norbital,nkamax+1,nkamax+1),dtype='float64')

#define uniform grid
grid=np.zeros((nkamax+1,nkamax+1,2))
B=np.array([[2*pi/(a),2*pi/(a)],[2*pi/(sqrt(3)*a),-2*pi/(sqrt(3)*a)]]) #Transformation matrix
for nk1 in range(0,nkamax+1):
    for nk2 in range(0,nkamax+1):
        [[grid[nk1,nk2,0]],[grid[nk1,nk2,1]]]=B@[[(nk1)/nkamax],[(nk2)/nkamax]]-np.array([[2*pi/(a)],[0]])
#plt.plot(grid[:,:,0],grid[:,:,1],'x')
#plt.show()
#check connection
for nk1 in range(0,nkamax+1):
    for nk2 in range(0,nkamax+1):
        k=grid[nk1,nk2,:]
        E[:,nk1,nk2],V[:,:,nk1,nk2]=LA.eigh(H(k[0],k[1]))
        for nu in range(0,Norbital):
            for mu in range(0,Norbital):
                p[nu,mu,nk1,nk2,0]=me/hbar*(np.transpose(np.conj(V[:,nu,nk1,nk2]))@Vpkx(k[0],k[1])@V[:,mu,nk1,nk2])
                p[nu,mu,nk1,nk2,1]=me/hbar*(np.transpose(np.conj(V[:,nu,nk1,nk2]))@Vpky(k[0],k[1])@V[:,mu,nk1,nk2])
                if np.abs(E[nu,nk1,nk2]-E[mu,nk1,nk2])<1e-2:
                    Delta[nu,mu]+=1
print("momentum p: defined")
for nu in range(0,Norbital):
    for mu in range(0,Norbital):
        if Delta[nu,mu]!=0:
            Delta[nu,mu]=1
print("Connect Matrix Delta: OK")
for nk1 in range(0,nkamax+1):
    for nk2 in range(0,nkamax+1):
        k=grid[nk1,nk2,:] #E already had
        for nu in range(0,Norbital):
            for mu in range(0,Norbital):
                if Delta[nu,mu]==0:
                    Xi[nu,mu,nk1,nk2,0]=-1j*hbar/me*(np.transpose(np.conj(V[:,nu,nk1,nk2]))@(Vpkx(k[0],k[1]))@V[:,mu,nk1,nk2])/(E[nu,nk1,nk2]-E[mu,nk1,nk2])
                    Xi[nu,mu,nk1,nk2,1]=-1j*hbar/me*(np.transpose(np.conj(V[:,nu,nk1,nk2]))@(Vpky(k[0],k[1]))@V[:,mu,nk1,nk2])/(E[nu,nk1,nk2]-E[mu,nk1,nk2])
                    BerryCurvz[nu,nk1,nk2]+=-2*np.imag((np.transpose(np.conj(V[:,nu,nk1,nk2]))@(Vpkx(k[0],k[1]))@V[:,mu,nk1,nk2])*(np.transpose(np.conj(V[:,mu,nk1,nk2]))@(Vpky(k[0],k[1]))@V[:,nu,nk1,nk2]))/(E[nu,nk1,nk2]-E[mu,nk1,nk2])**2
                    #note: Xi intra not coded yes, Maybe can use formular in thong's article
print("\u03BE_{inter}, \u03A9_{inter}$: calculated")
#@jit(forceobj=True)
def RHSsbeVG(R,nt):
    C=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex128')
    for nu in range(0,Norbital):
        for mu in range(0,Norbital):
            for nu1 in range(0,Norbital):
                C[nu,mu,:,:,0]+=np.multiply(p[nu,nu1,:,:,0],rho[nu1,mu,:,:])-np.multiply(rho[nu,nu1,:,:],p[nu1,mu,:,:,0])
                C[nu,mu,:,:,1]+=np.multiply(p[nu,nu1,:,:,1],rho[nu1,mu,:,:])-np.multiply(rho[nu,nu1,:,:],p[nu1,mu,:,:,1])
    for nk1 in range(0,nkamax+1):
        for nk2 in range(0,nkamax+1):
            for nu in range(0,Norbital):
                for mu in range(nu,Norbital): #calculate RHS for uppper triangle
                    if np.divmod(nt,int(nt))[1]==0: 
                        RHS[nu,mu,nk1,nk2]=-1j/hbar*(E[nu,nk1,nk2]-E[mu,nk1,nk2])*R[nu,mu,nk1,nk2]-1j*qe/(hbar*me)*(Alec[nt,:]@C[nu,mu,nk1,nk2,:])
                    else:
                        t=nt*deltat+tmin
                        RHS[nu,mu,nk1,nk2]=-1j/hbar*(E[nu,nk1,nk2]-E[mu,nk1,nk2])*R[nu,mu,nk1,nk2]-1j*qe/(hbar*me)*(Alec[math.floor((t-tmin)/deltat)-1,:]@C[nu,mu,nk1,nk2,:]-Elec[math.floor((t-tmin)/deltat),:]@C[nu,mu,nk1,nk2,:]*(t-math.floor((t-tmin)/deltat)*deltat))
                    RHS[mu,nu,nk1,nk2]=np.conj(RHS[nu,mu,nk1,nk2]) #calculate under triangle
    return RHS
#@jit(forceobj=True)
def RK(R,nt):
    RKdummy[:,:,:,:,0]=RHSsbeVG(R,nt)
    RKdummy[:,:,:,:,1]=RHSsbeVG(R+RKdummy[:,:,:,:,0]*deltat/3,nt+1/3)
    RKdummy[:,:,:,:,2]=RHSsbeVG(R-RKdummy[:,:,:,:,0]*deltat/3+RKdummy[:,:,:,:,1]*deltat,nt+2/3)
    RKdummy[:,:,:,:,3]=RHSsbeVG(R+RKdummy[:,:,:,:,0]*deltat-RKdummy[:,:,:,:,1]*deltat+RKdummy[:,:,:,:,2]*deltat,nt+1)
    R+=deltat/8*(RKdummy[:,:,:,:,0]+3*RKdummy[:,:,:,:,1]+3*RKdummy[:,:,:,:,2]+RKdummy[:,:,:,:,3])
    return R
#calculation:
for nt in range(0,int((tmax-tmin)/deltat)-1):
    start=time.time()
    if nt!=0:
        rho=RK(rho,nt)
    else:
        if Norbital==6:
            rho[3,3,:,:]=1
            rho[0,0,:,:]=1
        if Norbital==3:
            rho[0,0,:,:]==1
    for nu in range(0,Norbital):
        for nk1 in range(0,nkamax+1):
            for nk2 in range(0,nkamax+1):
                N[nu,nt]+=np.real(rho[nu,nu,nk1,nk2])/((nkamax+1)**2)
                J[nt,:,0]+=((p[nu,nu,nk1,nk2,:]-qe*Alec[nt,:])*rho[nu,nu,nk1,nk2])/me #intra current
                for mu in range(0,Norbital):
                    if mu != nu:
                        J[nt,:,1]+=((p[nu,mu,nk1,nk2,:]-qe*Alec[nt,:])*rho[nu,mu,nk1,nk2])/me #inter current
        #checkin section
        if N[nu,nt]<0:
            print("Density of a band can't be negative")
            exit()
        if N[nu,nt]>(nkamax+1)**2:
            print("Density of a band can't over 1")
            exit()
    J[nt,:,2]+=J[nt,:,1]+J[nt,:,0] #total current
    writercur = csv.writer(open('Currentdensity.csv', 'a', newline =''), quoting = csv.QUOTE_ALL,delimiter =';')
    writercur.writerows(np.transpose([[nt*deltat+tmin],[J[nt,0,0]],[J[nt,1,0]],[J[nt,0,1]], [J[nt,1,1]],[J[nt,0,2]],[J[nt,1,2]]]))
    if Norbital==3:
        writerden = csv.writer(open('Density.csv', 'a', newline =''), quoting = csv.QUOTE_ALL,delimiter =';')
        writerden.writerows(np.transpose([[nt*deltat+tmin],[rho[0,0,0,0]], [rho[1,1,0,0]],[rho[2,2,0,0]]]))
    if Norbital==6:
        writerden = csv.writer(open('Density.csv', 'a', newline =''), quoting = csv.QUOTE_ALL,delimiter =';')
        writerden.writerows(np.transpose([[nt*deltat+tmin],[rho[0,0,0,0]], [rho[1,1,0,0]],[rho[2,2,0,0]],[rho[3,3,0,0]], [rho[4,4,0,0]],[rho[5,5,0,0]]]))
    end=time.time()
    print("Calculationed:", "%.2f" % (nt*deltat+tmin),"on",tmax,".deltat for a step=:",(end-start),".Estimate ends:", time.ctime(end+(end-start)*(int((tmax-tmin)/deltat)-nt)))
#fig, axs = plt.subplots(2, sharex=True)
#fig.suptitle('Sharing both axes')
#axs[0].plot(Elec[:,0],label='E(t)')
#axs[0].plot(Alec[:,0],label='A(t)')
#axs[1].plot(J[:,0,0],'--b',label="Intra")
#axs[1].plot(J[:,0,1],'--',label="Inter")
#axs[1].plot(J[:,0,2],'-r',label="Total")
#axs[0].set_ylabel(r'Electric Field[V/nm]')
#axs[1].legend()
#axs[1].set_ylabel(r'Current Density')
#plt.show()
#Jf=np.zeros((int(NomegaL/deltaomega),2,3),dtype='complex128')
#Jffreq=np.zeros((int(NomegaL/deltaomega)),dtype='float32')
#for i in range(0,int(NomegaL/deltaomega)):
#    for nt in range(0,int((tmax-tmin)/deltat)):
#        Jf[i,:,0]+=J[nt,:,0]*exp(-1j*i*deltaomega*omegaL*(nt*deltat+tmin-tpeak))*deltat
#        Jf[i,:,1]+=J[nt,:,1]*exp(-1j*i*deltaomega*omegaL*(nt*deltat+tmin-tpeak))*deltat
#        Jf[i,:,2]+=J[nt,:,2]*exp(-1j*i*deltaomega*omegaL*(nt*deltat+tmin-tpeak))*deltat
#        Jffreq[i]=i*deltaomega
#    print(i,"on", int(NomegaL/deltaomega))
#plt.plot(Jffreq,Jf[:,0,0],'x')
#plt.plot(Jffreq,Jf[:,0,1],'b--')
#plt.show()
