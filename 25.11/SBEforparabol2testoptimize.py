import numpy as np
from scipy.linalg import svd
from numpy import linalg as LA
from numpy import sin, cos, sqrt, pi, exp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
tmin=-200 #fs
tmax=200 #fs
nkamax=999
E0=2.2e-1 #V/nm
TauL=50 #fs
hbar=0.6582 #[eV.fs]
epsilonL=1.6 #[eV]
omegaL=epsilonL/hbar #[1/fs]
qe=-1 #[unit charge]
me=9.1094/1.6022 #[eV.fs^2/nm^2];
deltat=0.02 #fs
mc=0.067*me
mv=-0.46*me
Eg=1.5 #[eV]
Norbital=2
kmin=-3 #[nm^(-1)]
kmax=3 #[nm^(-1)]
deltak=(kmax-kmin)/nkamax #[nm^(-1)]
#grid
grid=np.linspace(kmin,kmax,num=nkamax+1)
rho=np.zeros((Norbital,Norbital,nkamax+1),dtype='complex128')
grad=np.zeros((Norbital,Norbital,nkamax+1),dtype='complex128')
rhodummy=np.zeros((Norbital,Norbital,nkamax+1),dtype='complex128')
RKdummy=np.zeros((Norbital,Norbital,nkamax+1,4),dtype='complex128')
E=np.zeros((2,nkamax+1),dtype="float64")
Elec=np.zeros(int((tmax-tmin)/deltat),dtype="float64")
P=np.zeros((int((tmax-tmin)/deltat)),dtype='float64')
N=np.zeros((int((tmax-tmin)/deltat)),dtype='float64')
for nk in range(0,nkamax+1):
    k=grid[nk]
    E[0,nk]=(hbar*k)**2/(2*mv)
    E[1,nk]=(hbar*k)**2/(2*mc)+Eg
for nt in range(0,int((tmax-tmin)/deltat)):
    t=nt*deltat
    Elec[nt]=E0*exp(-(tmin+t)**2/TauL**2)*cos(omegaL*(tmin+t))
Xi=np.array([[0,0.3],[0.3,0]])
def RHS(Rho,t):
    C=np.zeros((Norbital,Norbital,nkamax+1),dtype='complex128')
    for nk in range(0,nkamax+1):
        grad[:,:,nk]=(-Rho[:,:,(nk+2)%nkamax]+8*(Rho[:,:,(nk+1)%nkamax])-8*Rho[:,:,(nk-1)%nkamax]+Rho[:,:,(nk+2)%nkamax])/(12*deltak)
    for nu in range(0,Norbital):
        for mu in range(0,Norbital):
            for nu1 in range(0,Norbital):
                C[nu,mu,:]+=Xi[nu,nu1]*Rho[nu1,mu,:]-Rho[nu,nu1,:]*Xi[nu1,mu]
    for nu in range(0,Norbital):
        for mu in range(0,Norbital):
            rhodummy[nu,mu,:]=-1j/hbar*(E[nu,:]-E[mu,:])*Rho[nu,mu,:]+(qe/hbar)*(grad[nu,mu,:]-1j*C[nu,mu,:])*E0*exp(-(tmin+t)**2/TauL**2)*cos(omegaL*(tmin+t))
    return rhodummy
def RK(Rho,t):
    #tinh k1
    RKdummy[:,:,:,0]=RHS(Rho,t) #Rkdummy: k1
    RKdummy[:,:,:,1]=RHS(Rho+RKdummy[:,:,:,0]*deltat/3,t+deltat/3)
    RKdummy[:,:,:,2]=RHS(Rho-RKdummy[:,:,:,0]*deltat/3+RKdummy[:,:,:,1]*deltat,t+2*deltat/3)
    RKdummy[:,:,:,3]=RHS(Rho+RKdummy[:,:,:,0]*deltat-RKdummy[:,:,:,1]*deltat+RKdummy[:,:,:,2]*deltat,t+deltat)
    return Rho+deltat/8*(RKdummy[:,:,:,0]+3*RKdummy[:,:,:,1]+3*RKdummy[:,:,:,2]+RKdummy[:,:,:,3])
rho[0,0,:]=1
for nt in range(0,int((tmax-tmin)/deltat)):
    t=deltat*nt
    start=time.time()
    rho=RK(rho,t)
    if np.divmod(nt*deltat,10)[1]==0 and (nt*deltat+tmin)>=-20:
        plt.plot(grid,rho[1,1,:],label='nt*deltat+tmin')
        plt.legend()
    for nk in range(0,nkamax+1):
        P[nt]+=Xi[1,0]*(rho[0,1,nk]+np.conj(rho[0,1,nk]))
        N[nt]+=rho[1,1,nk]
    if N[nt]<0:
        print("Density can't be negative")
        exit()
    print(N[nt],P[nt])
    end=time.time()
    print("Calculationed:", "%.2f" % (nt*deltat+tmin),"on",tmax,".deltat for a step=:",(end-start),".Estimate ends:", time.ctime(end+(end-start)*(int((tmax-tmin)/deltat)-nt)))
plt.show()
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
plt.subplot(311)
plt.plot(Elec)
plt.subplot(312)
plt.plot(P)
plt.subplot(313)
plt.plot(N)
plt.show()