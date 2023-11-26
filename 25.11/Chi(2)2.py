import sys
sys.path.append('Z:/root/TB')
from variablei import variable,nkamax,Norbital
from parameter import para
import numpy as np
from numpy import linalg as LA
from numpy import sin, cos, sqrt, pi, exp
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 
plt.rcParams.update({'font.size': 20})
from numba import njit
import time
#time in fs
phi=0 #rad
theta=0 #rad
TauL=60*sqrt(2) #fs (the TauL is not have the form x^2/2\tauL**2)
hbar=0.6582 #[eV.fs]
epsilonL=1.81 #[eV]
deltaomega=0.02
omegaL=epsilonL/(hbar)# [1/fs]
qe=1 #[unit charge]
me=9.1094/1.6022 #[eV.fs^2/nm^2]; E=mc^2;[kg]=[J(s/m)^2]=1/qe [eV (s/m)^2]
Ndimension=2
#dephasing
T2=5.2 #[fs]
#fourier transform for harmonic
p=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex64')
Xi=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex64')
V=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex64')
rho=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
RKdummy=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,4),dtype='complex128')
BerryCurvz=np.zeros((Norbital,nkamax+1,nkamax+1),dtype='float64')
Delta=np.zeros((Norbital,Norbital),dtype='int16')
deltaE=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='float32')
E=np.zeros((Norbital,nkamax+1,nkamax+1),dtype='float64')
Chi2dummy=np.zeros((6,Ndimension,Ndimension,Ndimension,int(4/deltaomega)),dtype='complex')
epsilon=np.zeros((int(4/deltaomega)),dtype='float')
fig, axs = plt.subplots(2,sharex=True,figsize=(9, 11))
for j in range(int(4/deltaomega)):
    epsilon[j]=j*deltaomega
print(epsilon)
for i in range(6):
    start=time.time()
    print(i)
    a=para(i)[0]
    dS=((4/sqrt(3))*pi/(a*(nkamax+1)))**2*cos(np.arctan(1/sqrt(3)))
    grid,p,Xi,BerryCurvz,Delta,deltaE,E,V=variable(i)
    print("Variable: Defined")
    if Norbital==3:
        valence=[0]
        conduct=[1,2]
    elif Norbital==6:
        valence=[0,3]
        conduct=[1,2,4,5]
    print("Valance,conduct index: defined")
    for nk1 in range(nkamax+1):
        for nk2 in range(nkamax+1):
            k=grid[nk1,nk2,:]
    def Chi2(omega):
        gamma=1e-1
        Chi2=np.zeros((Ndimension,Ndimension,Ndimension),dtype='complex')
        dummy2=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,Ndimension,Ndimension),dtype='complex')
        for m in range(Ndimension):
            for n in range(Ndimension):
                dummy2[...,m,n]+=-np.einsum('ij...,jk...->ik...',p[...,m]/(deltaE+1j*gamma),p[...,n])+np.einsum('ij...,jk...->ik...',p[...,n],p[...,m]/(deltaE+1j*gamma))
        for i in range(Ndimension):
            for j in range(Ndimension):
                for k in range(Ndimension):
                    for ci in conduct:
                        for vi in valence:
                            Chi2[i,j,k]+=(qe**3/(hbar*me)**2)*np.sum(np.multiply(p[ci,vi,:,:,k]/(deltaE[ci,vi,:,:]/hbar-omega-1j*gamma),dummy2[vi,ci,:,:,i,j]))*dS*1.6
        return Chi2
    for j in range(int(4/deltaomega)):
        Chi2dummy[i,:,:,:,j]=Chi2(epsilon[j]/hbar)
    end=time.time()
    print(i, "cost:", end-start)
axs[0].plot(epsilon,np.real(Chi2dummy[0,1,0,0,:]),label='MoS2')
axs[0].plot(epsilon,np.real(Chi2dummy[1,1,0,0,:]),label='WS2')
axs[0].plot(epsilon,np.real(Chi2dummy[2,1,0,0,:]),label='MoSe2')
axs[0].plot(epsilon,np.real(Chi2dummy[3,1,0,0,:]),label='WSe2')
axs[0].plot(epsilon,np.real(Chi2dummy[4,1,0,0,:]),label='MoTe2')
axs[0].plot(epsilon,np.real(Chi2dummy[5,1,0,0,:]),label='WTe2')
axs[0].legend(fontsize='x-small')
axs[0].set_ylabel(r'$Re\{\sigma^{(2)}_{yxx}\} \: 10^{17} [m. A/(Vs)^2]$ ',fontsize='x-small')
axs[0].grid()
axs[0].set_title("(a)", loc='left',size=20)

axs[1].plot(epsilon,np.imag(Chi2dummy[0,1,0,0,:]),label='MoS2')
axs[1].plot(epsilon,np.imag(Chi2dummy[1,1,0,0,:]),label='WS2')
axs[1].plot(epsilon,np.imag(Chi2dummy[2,1,0,0,:]),label='MoSe2')
axs[1].plot(epsilon,np.imag(Chi2dummy[3,1,0,0,:]),label='WSe2')
axs[1].plot(epsilon,np.imag(Chi2dummy[4,1,0,0,:]),label='MoTe2')
axs[1].plot(epsilon,np.imag(Chi2dummy[5,1,0,0,:]),label='WTe2')
axs[1].set_xlabel("Energy of Photon [eV]")
axs[1].set_ylabel(r'$Im\{\sigma^{(2)}_{yxx}\} \: 10^{17} [m. A/(Vs)^2]$ ',fontsize='x-small')
axs[1].grid()
axs[1].set_title("(b)", loc='left',size=20)


#enlarge picture
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
fig.savefig('Chi2new.svg', format='svg', dpi=1200)
fig.savefig('Chi2new.eps', format='eps')
#plt.show()