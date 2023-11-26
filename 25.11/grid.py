import numpy as np
import time
from scipy import constants as Const
from numpy import sin,cos,sqrt,pi
import matplotlib.pyplot as plt
from hamud import Norbital,a
from variableud import variable
plt.rcParams.update({'font.size': 20})
start=time.time()
nkamax=10
grid=np.zeros((nkamax+1,nkamax+1,2))
gridbigger=np.zeros((3*nkamax+3,3*nkamax+3,2))
B=np.array([[2*pi/(a),2*pi/(a)],[2*pi/(sqrt(3)*a),-2*pi/(sqrt(3)*a)]])

#B=np.array([[1,0],[0,1]])
for nk1 in range(0,nkamax+1):
    for nk2 in range(0,nkamax+1):
        [[grid[nk1,nk2,0]],[grid[nk1,nk2,1]]]=B@[[(nk1)/nkamax],[(nk2)/nkamax]]-np.array([[2*pi/(a)],[0]])
#plt.plot(grid[:,:,0],grid[:,:,1],'gx')
#for nk1 in range(0,3*nkamax+3):
#    for nk2 in range(0,3*nkamax+3):
#        [[gridbigger[nk1,nk2,0]],[gridbigger[nk1,nk2,1]]]=B@[[(nk1)/nkamax],[(nk2)/nkamax]]-np.array([[6*pi/(a)],[0]])
#plt.plot(gridbigger[(nkamax):2*(nkamax)+1,(nkamax):2*(nkamax)+1,0],gridbigger[(nkamax):2*(nkamax)+1,(nkamax):2*(nkamax)+1,1],'bx')
#plt.plot(gridbigger[0:(nkamax+1),0:(nkamax+1),0],grid[0:(nkamax+1),0:(nkamax+1),1],'mx')
#plt.plot(gridbigger[2*(nkamax):3*(nkamax)+1,2*(nkamax):3*(nkamax)+1,0],gridbigger[2*(nkamax):3*(nkamax)+1,2*(nkamax):3*(nkamax)+1,1],'mx')
#plt.plot(gridbigger[2*(nkamax):3*(nkamax)+1,0:nkamax+1,0],gridbigger[2*(nkamax):3*(nkamax)+1,0:nkamax+1,1],'rx')
#plt.plot(gridbigger[0:nkamax+1,2*(nkamax):3*(nkamax)+1,0],gridbigger[0:nkamax+1,2*(nkamax):3*(nkamax)+1,1],'rx')
#plt.plot(gridbigger[1*(nkamax):2*(nkamax)+1,2*(nkamax):3*(nkamax)+1,0],gridbigger[1*(nkamax):2*(nkamax)+1,2*(nkamax):3*(nkamax)+1,1],'kx')
#plt.plot(gridbigger[0*(nkamax):1*(nkamax)+1,1*(nkamax):2*(nkamax)+1,0],gridbigger[0*(nkamax):1*(nkamax)+1,1*(nkamax):2*(nkamax)+1,1],'gx')
#plt.plot(gridbigger[1*(nkamax):2*(nkamax)+1,0*(nkamax):1*(nkamax)+1,0],gridbigger[1*(nkamax):2*(nkamax)+1,0*(nkamax):1*(nkamax)+1,1],'kx')
#plt.plot(gridbigger[2*(nkamax):3*(nkamax)+1,1*(nkamax):2*(nkamax)+1,0],gridbigger[2*(nkamax):3*(nkamax)+1,1*(nkamax):2*(nkamax)+1,1],'gx')
#high sysmestry point
#K
fig, axs = plt.subplots(2,figsize=(9, 11))
axs[0].plot(4*pi/(3*a),0,'ro',label='K')
#axs[0].annotate(r"K", # this is the text
#                 (4*pi/(3*a),0), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(5,0), # distance from text to points (x,y)
#                 ha='left')
axs[0].plot(-2*pi/(3*a),-2*pi/(sqrt(3)*a),'ro',)
#axs[0].annotate(r"K", # this is the text
#                 (-2*pi/(3*a),-2*pi/(sqrt(3)*a)), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(-5,-10), # distance from text to points (x,y)
#                 ha='right')    
axs[0].plot(-2*pi/(3*a),2*pi/(sqrt(3)*a),'ro',)
#axs[0].annotate(r"K", # this is the text
#                 (-2*pi/(3*a),2*pi/(sqrt(3)*a)), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(0,5), # distance from text to points (x,y)
#                 ha='left')
#-K
axs[0].plot(-4*pi/(3*a),0,'bo',label="K'")
#axs[0].annotate(r"K'", # this is the text
#                 (-4*pi/(3*a),0), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(10,10), # distance from text to points (x,y)
#                 ha='left')
axs[0].plot(2*pi/(3*a),-2*pi/(sqrt(3)*a),'bo')
#axs[0].annotate(r"K'", # this is the text
#                 (2*pi/(3*a),-2*pi/(sqrt(3)*a)), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(3,-10), # distance from text to points (x,y)
#                 ha='left')
axs[0].plot(2*pi/(3*a),2*pi/(sqrt(3)*a),'bo')
#axs[0].annotate(r"K'", # this is the text
#                 (2*pi/(3*a),2*pi/(sqrt(3)*a)), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(0,5), # distance from text to points (x,y)
#                 ha='left')
#M
axs[0].plot(pi/(a),pi/(sqrt(3)*a),'gs',label='M')
#axs[0].annotate(r"M", # this is the text
#                 (pi/a,pi/(sqrt(3)*a)), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(0,10), # distance from text to points (x,y)
#                 ha='left')
axs[0].plot(2*pi/(a),0,'gs')
#axs[0].annotate(r"M", # this is the text
#                 (2*pi/a,0), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(0,10), # distance from text to points (x,y)
#                 ha='left')
axs[0].plot(-2*pi/(a),0,'gs')
#axs[0].annotate(r"M", # this is the text
#                 (-2*pi/a,0), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(0,10), # distance from text to points (x,y)
#                 ha='center')#
axs[0].plot(0,-2*pi/(sqrt(3)*a),'gs')
#axs[0].annotate(r"M", # this is the text
#                 (0,-2*pi/(sqrt(3)*a)), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(0,-20), # distance from text to points (x,y)
#                 ha='right')
axs[0].plot(0,2*pi/(sqrt(3)*a),'gs')
#axs[0].annotate(r"M", # this is the text
#                 (0,2*pi/(sqrt(3)*a)), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(0,5), # distance from text to points (x,y)
#                 ha='right')
#line
axs[0].plot([4*pi/(3*a),2*pi/(3*a)],[0,2*pi/(sqrt(3)*a)],'g-')
axs[0].plot([-2*pi/(3*a),2*pi/(3*a)],[2*pi/(sqrt(3)*a),2*pi/(sqrt(3)*a)],'g-')
axs[0].plot([-2*pi/(3*a),-4*pi/(3*a)],[2*pi/(sqrt(3)*a),0],'g-')
axs[0].plot([-2*pi/(3*a),-4*pi/(3*a)],[-2*pi/(sqrt(3)*a),0],'g-')
axs[0].plot([-2*pi/(3*a),2*pi/(3*a)],[-2*pi/(sqrt(3)*a),-2*pi/(sqrt(3)*a)],'g-')
axs[0].plot([4*pi/(3*a),2*pi/(3*a)],[0,-2*pi/(sqrt(3)*a)],'g-')
axs[0].plot([2*pi/(a),0],[0,2*pi/(sqrt(3)*a)],'r--')
axs[0].plot([-2*pi/(a),0],[0,2*pi/(sqrt(3)*a)],'r--')
axs[0].plot([2*pi/(a),0],[0,-2*pi/(sqrt(3)*a)],'r--')
axs[0].plot([-2*pi/(a),0],[0,-2*pi/(sqrt(3)*a)],'r--')
axs[0].arrow(-2*pi/(a),0, 2*pi/(a)/10,2*pi/(sqrt(3)*a)/10,width=0.5)
#axs[0].annotate(r"$u_1$", # this is the text
#                 (-2*pi/(a)*0.9,2*pi/(a)/10*1.1), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(10,-5), # distance from text to points (x,y)
#                 ha='center')
#axs[0].annotate(r"$u_2$", # this is the text
#                 (-2*pi/(a)*0.9,-2*pi/(a)/10*1.1), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(10,-5), # distance from text to points (x,y)
#                 ha='center')
axs[0].arrow(-2*pi/(a),0, 2*pi/(a)/10,-2*pi/(sqrt(3)*a)/10,width=0.5)
#axs[0].annotate(r"$k_x$", # this is the text
#                 (2,0), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(10,-10), # distance from text to points (x,y)
#                 ha='left')
#axs[0].annotate(r"$k_y$", # this is the text
#                 (0,2), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(-9,5), # distance from text to points (x,y)
#                 ha='center')
#axs[0].annotate("O", # this is the text
#                 (0,0), # these are the coordinates to position the label
#                 textcoords="offset points", # how to position the text
#                 xytext=(-10,-6), # distance from text to points (x,y)
#                 ha='right')
axs[0].arrow(0,0,2,0, width=0.5)
axs[0].arrow(0,0,0,2, width=0.5)
#xticks
#axs[0].xticks([-6*pi/(a),-2*pi/(a),-4*pi/(3*a), 0,4*pi/(3*a), 2*pi/(a),6*pi/(a)], [r'$-\frac{6\pi}{a}$',r'$-\frac{2\pi}{a}$',r'$-\frac{4\pi}{3a}$','0',r'$\frac{4\pi}{3a}$',r'$\frac{2\pi}{a}$',r'$\frac{6\pi}{a}$'])
#axs[0].set_xticks([-6*pi/(a),-2*pi/(a), 0, 2*pi/(a),6*pi/(a)], [r'$-\frac{6\pi}{a}$',r'$-\frac{2\pi}{a}$','0',r'$\frac{2\pi}{a}$',r'$\frac{6\pi}{a}$'])
#axs[0].set_yticks([-6*pi/(sqrt(3)*a),-2*pi/(sqrt(3)*a), 0, 2*pi/(sqrt(3)*a),6*pi/(sqrt(3)*a)], [r'$-\frac{6\pi}{\sqrt{3}a}$',r'$-\frac{2\pi}{\sqrt{3}a}$', 0, r'$\frac{2\pi}{\sqrt{3}a}$',r'$\frac{6\pi}{\sqrt{3}a}$'])
axs[0].set_yticks([-2*pi/(sqrt(3)*a), 0, 2*pi/(sqrt(3)*a)], [r'$-\frac{2\pi}{\sqrt{3}a}$', 0, r'$\frac{2\pi}{\sqrt{3}a}$'])
#plt.title("Multies Unit Cells In Reciprocal Space", y=-0.1)
axs[0].plot(grid[:,:,0],grid[:,:,1],'b')
axs[0].plot(grid[:,:,0],np.transpose(grid[:,:,1],(1,0)),'b')
axs[0].grid()
axs[0].legend(fontsize='x-small')
axs[0].set_title("(a)", loc='left',size=15)
axs[0].set_xticks([],[])
from variableud import nkamax
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
print("Variable: Defined")
dummy1=np.zeros((Norbital,nkamax+1),dtype='complex128')
dummy2=np.zeros((nkamax+1,2),dtype='complex128')
#RHS=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
for nk in range(nkamax+1):
    dummy2[nk,:]=grid[nk,nk,:]
    for nu in range(Norbital):
        dummy1[nu,nk]=E[nu,nk,nk]
for i in range(Norbital):
    axs[1].plot(dummy2[:,0],dummy1[i,:],'r-')
#axs[1].set_xticks([-6*pi/(a),-2*pi/(a), 0, 2*pi/(a),6*pi/(a)], [r'$-\frac{6\pi}{a}$',r'$-\frac{2\pi}{a}$','0',r'$\frac{2\pi}{a}$',r'$\frac{6\pi}{a}$'])
axs[1].set_xticks([dummy2[0,0],dummy2[int((nkamax)/6),0],dummy2[int((nkamax)/2),0],dummy2[int(5*(nkamax)/6),0],dummy2[int(nkamax),0]],["M","-K",r"$\Gamma$","K","M"])
axs[1].set_title("(b)", loc='left',size=15)
plt.grid()
plt.ylabel("Energy [eV]")
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
fig.savefig('GridandBS.svg', format='svg', dpi=1200)
fig.savefig('GridandBS.eps', format='eps')
#plt.show()
#end=time.time()
#print(end-start)
#plt.show()