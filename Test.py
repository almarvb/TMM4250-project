import numpy as np


num_dofs = 27
uVec   = np.arange(1,27)
#uVec = uVec.T
print ( uVec)
for iel in range(8):

    uLiten =uVec[3*iel:3*iel+6]
    print (uLiten)