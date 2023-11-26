import sys
sys.path.append('Z:/root/TB')
import matplotlib
matplotlib.use('TkAgg')
from variableud import variable,nkamax,Norbital, dS
import numpy as np
from numpy import linalg as LA
from numpy import sin, cos, sqrt, pi, exp
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
from numba import jit
#time in fs
phi=0 #rad
theta=0 #rad
TauL=60*sqrt(2) #fs (the TauL is not have the form x^2/2\tauL**2)
hbar=0.6582 #[eV.fs]
epsilonL=1.81 #[eV]
deltaomega=0.01
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
grid,p,Xi,BerryCurvz,Delta,deltaE,E,V=variable()
epsilon=np.zeros((int(4/deltaomega)),dtype='float')
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
    for a in range(Ndimension):
        for b in range(Ndimension):
            for ci in conduct:
                for vi in valence:
                   for la in range(Norbital):
                        if la !=ci:
                            dummy2[vi,ci,:,:,a,b]+=-p[vi,la,:,:,b]*p[la,ci,:,:,a]/(deltaE[la,ci,:,:]/hbar+1j*gamma)
                        if la !=vi:
                            dummy2[vi,ci,:,:,a,b]+=p[vi,la,:,:,a]/(deltaE[vi,la,:,:]/hbar+1j*gamma)*p[la,ci,:,:,b]
    for i in range(Ndimension):
        for j in range(Ndimension):
            for k in range(Ndimension):
                for ci in conduct:
                    for vi in valence:
                        Chi2[i,j,k]+=(qe**3/(hbar**2*me**3))*np.sum(np.multiply(p[ci,vi,:,:,k]/(deltaE[ci,vi,:,:]/hbar-omega-1j*gamma),dummy2[vi,ci,:,:,i,j]))*dS*1.6
    return Chi2
for j in range(int(4/deltaomega)):
    epsilon[j]=j*deltaomega
Chi2dummy=np.zeros((Ndimension,Ndimension,Ndimension,int(4/deltaomega)),dtype='complex')
for j in range(int(4/deltaomega)):
    print(epsilon[j])
    Chi2dummy[...,j]=Chi2(epsilon[j]/hbar)
fig, axs = plt.subplots(2,sharex=True,figsize=(9, 11))
axs[0].plot(epsilon,np.real(Chi2dummy[0,0,1,:]),'b',label=r'$\sigma_{xxy}$')
axs[0].plot(epsilon,np.real(Chi2dummy[0,1,0,:]),'y',label=r'$\sigma_{xyx}$')
axs[0].plot(epsilon,np.real(Chi2dummy[1,0,0,:]),'r',label=r'$\sigma_{yxx}$')
axs[0].plot(epsilon,np.real(Chi2dummy[1,1,1,:]),'g',label=r'$\sigma_{yyy}$')
axs[0].set_xlabel("Energy of Photon [eV]")
axs[0].set_ylabel(r'$Re\{\sigma^{(2)}_{yxx}\} \: 10^{17} [m. A/(Vs)^2]$ ')
axs[0].legend(fontsize='x-small')
axs[0].grid()
axs[1].plot(epsilon,np.imag(Chi2dummy[0,0,1,:]),'b',label=r'$\sigma_{xxy}$')
axs[1].plot(epsilon,np.imag(Chi2dummy[0,1,0,:]),'y',label=r'$\sigma_{xyx}$')
axs[1].plot(epsilon,np.imag(Chi2dummy[1,0,0,:]),'r',label=r'$\sigma_{yxx}$')
axs[1].plot(epsilon,np.imag(Chi2dummy[1,1,1,:]),'g',label=r'$\sigma_{yyy}$')
axs[1].set_ylabel(r'$Im\{\sigma^{(2)}_{yxx}\} \: 10^{17} [m. A/(Vs)^2]$ ')
axs[1].legend(fontsize='x-small')
axs[1].grid()
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
plt.savefig('images/Chi2MoS2full.svg', format='svg', dpi=1200)
plt.savefig('images/Chi2MoS2full.eps', format='eps')