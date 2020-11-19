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
'''
num_nodes = 9
beamModel = sbeam.CantileverWithEndMoment(num_nodes)
#beamModel = CantileverWithEndMoment(num_nodes)

sbeam.solveLinearSteps(beamModel)
#sbeam.solveArchLength(beamModel)

num_steps = len(beamModel.load_history)

for iStep in range(num_steps):
   print("LoadFactor= {:12.3e}".format(beamModel.load_history[iStep]))
   print("dispVec={:}".format(iStep))
   print(beamModel.disp_history[iStep])

step_inc = (num_steps // 10)
for iStep in range(0,len(beamModel.load_history), step_inc):
    beamModel.plotDispState(iStep)

print("End")
'''

#--------------------------------------------------------------------------------------------
# ------------------ Preform non-linear solution (load controll + Newton itterations )
num_nodes = 9
beamModel = sbeam.CantileverWithEndMoment(num_nodes)
load_steps=0.01 
N_steps=200 
max_iter=30
#beamModel = CantileverWithEndMoment(num_nodes)

sbeam.solveNonlinLoadControl(beamModel,load_steps, N_steps, max_iter)
#Nonlinear load Control, with newton steps

num_steps = len(beamModel.load_history)

for iStep in range(num_steps):
   print("LoadFactor= {:12.3e}".format(beamModel.load_history[iStep]))
   print("dispVec={:}".format(iStep))
   print(beamModel.disp_history[iStep])

step_inc = (num_steps // 10)
if step_inc < 1:
    step_inc = 1
for iStep in range(0,len(beamModel.load_history), step_inc):
    beamModel.plotDispState(iStep)

print("End")
