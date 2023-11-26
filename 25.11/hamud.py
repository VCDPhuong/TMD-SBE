import numpy as np
from numpy import sin, cos, sqrt, pi, exp
from numba import njit
from parameter import para
a,e1,e2,t0,t1,t2,t11,t12,t22,r0,r1,r2,r11,r12,u0,u1,u2,u11,u12,u22,lam=para(0)
#0 MoS2, 1 WS2, 2 MoSe2, 3 WSe2, 4 MoTe2, 5 WTe2
@njit
def H(kx,ky):
    h0=e1+2*t0*(2*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*kx*a/2))+2*r0*(2*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*ky*a*sqrt(3)/2))+2*u0*(2*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+cos(4*kx*a/2))
    h1=complex(-2*sqrt(3.)*t2*sin(kx*a/2)*sin(ky*a*sqrt(3)/2)+2*(r1+r2)*sin(3*kx*a/2)*sin(ky*a*sqrt(3)/2)-2*sqrt(3.)*u2*sin(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2) , 2*t1*sin(kx*a/2)*(2*cos(kx*a/2)+cos(ky*a*sqrt(3)/2))+2*(r1-r2)*sin(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+2*u1*sin(2*kx*a/2)*(2*cos(2*kx*a/2)+cos(2*ky*a*sqrt(3)/2)))
    h2=complex(2*t2*(cos(2*kx*a/2)-cos(kx*a/2)*cos(ky*a*sqrt(3)/2))-2/sqrt(3.)*(r1 + r2)*(cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)-cos(2*ky*a*sqrt(3)/2))+2*u2*(cos(4*kx*a/2)-cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)),2*sqrt(3.)*t1*cos(kx*a/2)*sin(ky*a*sqrt(3)/2)+2/sqrt(3.)*sin(ky*a*sqrt(3)/2)*(r1-r2)*(cos(3*kx*a/2)+2*cos(ky*a*sqrt(3)/2))+2*sqrt(3.)*u1*cos(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2))
    h11=e2+(t11+3*t22)*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+2*t11*cos(2*kx*a/2)+4*r11*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+2*(r11+sqrt(3.)*r12)*cos(2*ky*a*sqrt(3)/2)+(u11+3*u22)*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+2*u11*cos(4*kx*a/2)
    h12=complex(sqrt(3.)*(t22-t11)*sin(kx*a/2.)*sin(ky*a*sqrt(3)/2.)+4*r12*sin(3*kx*a/2)*sin(ky*a*sqrt(3)/2)+sqrt(3.)*(u22-u11)*sin(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2),4*t12*sin(kx*a/2)*(cos(kx*a/2)-cos(ky*a*sqrt(3)/2))+4*u12*sin(2*kx*a/2)*(cos(2*kx*a/2)-cos(2*ky*a*sqrt(3)/2)))
    h22=e2+(3*t11+t22)*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+2*t22*cos(2*kx*a/2)+2*r11*(2*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*ky*a*sqrt(3)/2))+2/sqrt(3.)*r12*(4*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)-cos(2*ky*a*sqrt(3)/2))+(3*u11+u22)*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+2*u22*cos(4*kx*a/2)
    return [[h0,h1,h2,0,0,0],[np.conj(h1),h11,h12+1j*lam,0,0,0],[np.conj(h2),np.conj(h12)-1j*lam,h22,0,0,0],[0,0,0,h0,h1,h2],[0,0,0,np.conj(h1),h11,h12-1j*lam],[0,0,0,np.conj(h2),np.conj(h12)+1j*lam,h22]]
Norbital,col=np.shape(H(0,0))
@njit
def Hamu(kx,ky):
    h0=e1+2*t0*(2*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*kx*a/2))+2*r0*(2*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*ky*a*sqrt(3)/2))+2*u0*(2*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+cos(4*kx*a/2))
    h1=complex(-2*sqrt(3.)*t2*sin(kx*a/2)*sin(ky*a*sqrt(3)/2)+2*(r1+r2)*sin(3*kx*a/2)*sin(ky*a*sqrt(3)/2)-2*sqrt(3.)*u2*sin(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2) , 2*t1*sin(kx*a/2)*(2*cos(kx*a/2)+cos(ky*a*sqrt(3)/2))+2*(r1-r2)*sin(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+2*u1*sin(2*kx*a/2)*(2*cos(2*kx*a/2)+cos(2*ky*a*sqrt(3)/2)))
    h2=complex(2*t2*(cos(2*kx*a/2)-cos(kx*a/2)*cos(ky*a*sqrt(3)/2))-2/sqrt(3.)*(r1 + r2)*(cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)-cos(2*ky*a*sqrt(3)/2))+2*u2*(cos(4*kx*a/2)-cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)),2*sqrt(3.)*t1*cos(kx*a/2)*sin(ky*a*sqrt(3)/2)+2/sqrt(3.)*sin(ky*a*sqrt(3)/2)*(r1-r2)*(cos(3*kx*a/2)+2*cos(ky*a*sqrt(3)/2))+2*sqrt(3.)*u1*cos(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2))
    h11=e2+(t11+3*t22)*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+2*t11*cos(2*kx*a/2)+4*r11*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+2*(r11+sqrt(3.)*r12)*cos(2*ky*a*sqrt(3)/2)+(u11+3*u22)*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+2*u11*cos(4*kx*a/2)
    h12=complex(sqrt(3.)*(t22-t11)*sin(kx*a/2.)*sin(ky*a*sqrt(3)/2.)+4*r12*sin(3*kx*a/2)*sin(ky*a*sqrt(3)/2)+sqrt(3.)*(u22-u11)*sin(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2),4*t12*sin(kx*a/2)*(cos(kx*a/2)-cos(ky*a*sqrt(3)/2))+4*u12*sin(2*kx*a/2)*(cos(2*kx*a/2)-cos(2*ky*a*sqrt(3)/2)))
    h22=e2+(3*t11+t22)*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+2*t22*cos(2*kx*a/2)+2*r11*(2*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*ky*a*sqrt(3)/2))+2/sqrt(3.)*r12*(4*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)-cos(2*ky*a*sqrt(3)/2))+(3*u11+u22)*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+2*u22*cos(4*kx*a/2)
    return [[h0,h1,h2],[np.conj(h1),h11,h12+1j*lam],[np.conj(h2),np.conj(h12)-1j*lam,h22]]
@njit
def Hamd(kx,ky):
    h0=e1+2*t0*(2*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*kx*a/2))+2*r0*(2*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*ky*a*sqrt(3)/2))+2*u0*(2*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+cos(4*kx*a/2))
    h1=complex(-2*sqrt(3.)*t2*sin(kx*a/2)*sin(ky*a*sqrt(3)/2)+2*(r1+r2)*sin(3*kx*a/2)*sin(ky*a*sqrt(3)/2)-2*sqrt(3.)*u2*sin(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2) , 2*t1*sin(kx*a/2)*(2*cos(kx*a/2)+cos(ky*a*sqrt(3)/2))+2*(r1-r2)*sin(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+2*u1*sin(2*kx*a/2)*(2*cos(2*kx*a/2)+cos(2*ky*a*sqrt(3)/2)))
    h2=complex(2*t2*(cos(2*kx*a/2)-cos(kx*a/2)*cos(ky*a*sqrt(3)/2))-2/sqrt(3.)*(r1 + r2)*(cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)-cos(2*ky*a*sqrt(3)/2))+2*u2*(cos(4*kx*a/2)-cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)),2*sqrt(3.)*t1*cos(kx*a/2)*sin(ky*a*sqrt(3)/2)+2/sqrt(3.)*sin(ky*a*sqrt(3)/2)*(r1-r2)*(cos(3*kx*a/2)+2*cos(ky*a*sqrt(3)/2))+2*sqrt(3.)*u1*cos(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2))
    h11=e2+(t11+3*t22)*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+2*t11*cos(2*kx*a/2)+4*r11*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+2*(r11+sqrt(3.)*r12)*cos(2*ky*a*sqrt(3)/2)+(u11+3*u22)*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+2*u11*cos(4*kx*a/2)
    h12=complex(sqrt(3.)*(t22-t11)*sin(kx*a/2.)*sin(ky*a*sqrt(3)/2.)+4*r12*sin(3*kx*a/2)*sin(ky*a*sqrt(3)/2)+sqrt(3.)*(u22-u11)*sin(2*kx*a/2)*sin(2*ky*a*sqrt(3)/2),4*t12*sin(kx*a/2)*(cos(kx*a/2)-cos(ky*a*sqrt(3)/2))+4*u12*sin(2*kx*a/2)*(cos(2*kx*a/2)-cos(2*ky*a*sqrt(3)/2)))
    h22=e2+(3*t11+t22)*cos(kx*a/2)*cos(ky*a*sqrt(3)/2)+2*t22*cos(2*kx*a/2)+2*r11*(2*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)+cos(2*ky*a*sqrt(3)/2))+2/sqrt(3.)*r12*(4*cos(3*kx*a/2)*cos(ky*a*sqrt(3)/2)-cos(2*ky*a*sqrt(3)/2))+(3*u11+u22)*cos(2*kx*a/2)*cos(2*ky*a*sqrt(3)/2)+2*u22*cos(4*kx*a/2)
    return [[h0,h1,h2],[np.conj(h1),h11,h12-1j*lam],[np.conj(h2),np.conj(h12)+1j*lam,h22]]
@njit
def Vpkx(kx,ky):
    delta=1e-5
    return (np.array(H(kx+delta,ky), dtype='complex128')-np.array(H(kx-delta,ky),dtype='complex128'))/(2*delta)
@njit
def Vpky(kx,ky):
    delta=1e-5
    return (np.array(H(kx,ky+delta), dtype='complex128')-np.array(H(kx,ky-delta),dtype='complex128'))/(2*delta)