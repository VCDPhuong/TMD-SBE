import numpy as np  
from variablei import variable, nkamax
Norbital=6
E=np.zeros((Norbital,nkamax+1,nkamax+1),dtype='float64')
p=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex128')
Xi=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1,2),dtype='complex128')
Delta=np.zeros((Norbital,Norbital),dtype='int16')
grid=np.zeros((nkamax+1,nkamax+1,2))
V=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='complex128')
BerryCurvz=np.zeros((Norbital,nkamax+1,nkamax+1),dtype='float64')
deltaE=np.zeros((Norbital,Norbital,nkamax+1,nkamax+1),dtype='float32')
Gap=[6,6,6,6,6,6]
for i in range(6):
    grid,p,Xi,BerryCurvz,Delta,deltaE,E,V=variable(i)
    print(i)
    for nk1 in range(nkamax+1):
        for nk2 in range(nkamax+1):
            if (Gap[i]>deltaE[1,0,nk1,nk2]):
                Gap[i]=deltaE[1,0,nk1,nk2]
            elif (Gap[i]>deltaE[4,3,nk1,nk2]):
                Gap[i]=deltaE[4,3,nk1,nk2]
print(Gap)