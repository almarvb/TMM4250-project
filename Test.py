import numpy as np
import math

R = np.array([1,-2,0,0,0,0])

result = R.dot(R) < 1 

print (result)