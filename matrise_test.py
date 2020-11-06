import numpy as np
import math

disp_global = np.array([1,2,3,4,5,6])


diffx = disp_global[3]-disp_global[0]
diffy = disp_global[4]-disp_global[1]
Ld = math.sqrt(diffy**2+diffx**2)
#print(diffx)
#print(diffy)
print('Length_deformed is equal to: ' + str(Ld))