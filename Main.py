
import BeamModels as sbeam


'''This main script will solve the selected beam model using the selected solver and parameters.'''

'''Decide number of nodes in structure'''
num_nodes = 9

'''Select 1 beam model'''
beamModel = sbeam.CantileverWithEndMoment(num_nodes)
#beamModel = sbeam.SimplySupportedBeamModel(num_nodes)
#beamModel = sbeam.DeepArchModel(num_nodes)

'''Select applicable solver parameters'''
load_steps=0.01 #in case of LoadControl
arcLength=0.02 #in case of ArcLength
max_steps=100
max_iter=30

'''Select 1 solver method:'''
sbeam.solveNonlinLoadControl(beamModel,load_steps, max_steps, max_iter)
#sbeam.solveArcLength(beamModel, arcLength, max_steps, max_iter)
#sbeam.solveLinearSteps(beamModel, load_steps, max_steps)



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
#--------------------------------------------------------------------------------------------