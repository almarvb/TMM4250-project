#Cantilever model.
#-------------------------------------------------------------------------------------------
# Tenkte vi kunne bruke denne som et master script for å kjøre Cantiliver modellen vår.
# Da prøver vi å ikke endra sååå mye fra CorotBeam og TwoSimpleBeamModels men istede importere
# fra dem som pakker.
#-------------------------------------------------------------------------------------------
# Pakkene vi trenger
import math
import numpy as np
import matplotlib.pyplot as plt
import CorotBeam_with_TODO as CorotBeam
import matplotlib.animation as anm
from copy import deepcopy
import TwoSimpleBeamModels_with_TODO as sbeam

#--------------------------------------------------------------------------------------------
# ------------------ Perform linear solution

num_nodes = 9
beamModel = sbeam.CantileverWithEndMoment(num_nodes)
#beamModel = CantileverWithEndMoment(num_nodes)

sbeam.solveLinearSteps(beamModel)

num_steps = len(beamModel.load_history)

for iStep in range(num_steps):
   print("LoadFactor= {:12.3e}".format(beamModel.load_history[iStep]))
   print("dispVec={:}".format(iStep))
   print(beamModel.disp_history[iStep])

step_inc = (num_steps // 10)
for iStep in range(0,len(beamModel.load_history), step_inc):
    beamModel.plotDispState(iStep)

print("End")

#Eg skriv litt greier her for å sjå om det funker