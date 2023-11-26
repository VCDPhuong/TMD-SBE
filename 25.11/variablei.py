import numpy as np
from hamudi import Hamu,Hamd,H,Vpkx,Vpky,Norbital
from numpy import sin, cos, sqrt, pi, exp
from numpy import linalg as LA
from parameter import para
hbar=0.6582 #[eV.fs]
qe=1 #[unit charge]
me=9.1094/1.6022 #[eV.fs^2/nm^2]
nkamax=199
E=np.zeros((Norbital,nkamax+1,nkamax+1),dtype='float64')
p=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex128')
Xi=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex128')
Delta=np.zeros((Norbital,Norbital),dtype='int16')
grid=np.zeros((nkamax+1,nkamax+1,2))
V=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
BerryCurvz=np.zeros((Norbital,nkamax+1,nkamax+1),dtype='float64')
deltaE=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='float32')

def variable(i):
    a=para(i)[0]
    B=np.array([[2*pi/(a),2*pi/(a)],[2*pi/(sqrt(3)*a),-2*pi/(sqrt(3)*a)]]) #Transformation matrix
    for nk1 in range(0,nkamax+1):
        for nk2 in range(0,nkamax+1):
            [[grid[nk1,nk2,0]],[grid[nk1,nk2,1]]]=B@[[(nk1)/nkamax],[(nk2)/nkamax]]-np.array([[2*pi/(a)],[0]])
    for nk1 in range(nkamax+1):
        for nk2 in range(nkamax+1):
            k=grid[nk1,nk2]
            E[0:3,nk1,nk2],V[0:3,0:3,nk1,nk2]=LA.eigh(Hamu(i,k[0],k[1]))
            E[3:6,nk1,nk2],V[3:6,3:6,nk1,nk2]=LA.eigh(Hamd(i,k[0],k[1]))
    for nu in range(Norbital):
        for mu in range(Norbital):
            deltaE[nu,mu,:,:]=E[nu,:,:]-E[mu,:,:]
    for nk1 in range(nkamax+1):
        for nk2 in range(nkamax+1):
            k=grid[nk1,nk2,:]
            for nu in range(0,Norbital):
                for mu in range(0,Norbital):
                    p[nu,mu,nk1,nk2,0]=me/hbar*(np.transpose(np.conj(V[:,nu,nk1,nk2]))@Vpkx(i,k[0],k[1])@V[:,mu,nk1,nk2]) #[fs/nm*eV]
                    p[nu,mu,nk1,nk2,1]=me/hbar*(np.transpose(np.conj(V[:,nu,nk1,nk2]))@Vpky(i,k[0],k[1])@V[:,mu,nk1,nk2]) #[fs/nm*eV]
                    if np.abs(deltaE[nu,mu,nk1,nk2])<1e-2:
                        Delta[nu,mu]+=1
    for nu in range(0,Norbital):
        for mu in range(0,Norbital):
            if Delta[nu,mu]!=0:
                Delta[nu,mu]=1
#    Delta=np.identity(Norbital)

    for nu in range(0,Norbital):
        for mu in range(0,Norbital):
            if Delta[nu,mu]==0:
                Xi[nu,mu,:,:,0]=-1j*hbar/me*p[nu,mu,:,:,0]/deltaE[nu,mu,:,:]
                Xi[nu,mu,:,:,1]=-1j*hbar/me*p[nu,mu,:,:,1]/deltaE[nu,mu,:,:]
                BerryCurvz[nu,:,:]+=-2*(hbar/me)**2*np.imag(p[nu,mu,:,:,0]*p[nu,mu,:,:,1])/(deltaE[nu,mu,:,:])**2
                #note: Xi intra not coded yes, Maybe can use formular in thong's article
    return grid,p,Xi,BerryCurvz,Delta,deltaE,E,V