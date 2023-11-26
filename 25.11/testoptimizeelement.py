import numpy as np
qe=1
a = np.random.random((2,2,2,2,2))
Jdummy = np.zeros((2,2,2,2,2))
Jdummy2 = np.zeros((2,2,2,2,2))
Alec=[[1,2]]
Jdummy=a[...,:]-qe*Alec[:]
for nu in range(0,2):
    for mu in range(0,2):
        for nk1 in range(2):
            for nk2 in range(2):
                Jdummy2[nu,mu,nk1,nk2,:]=a[nu,mu,nk1,nk2,:]-qe*Alec[:]
print(np.allclose(Jdummy,Jdummy2))