import math
import numpy as np
import matplotlib.pyplot as plt
import CorotBeam_with_TODO as CorotBeam
import matplotlib.animation as anm
from copy import deepcopy
import TwoSimpleBeamModels_with_TODO as sbeam


class DeepArchModel(sbeam.BeamModel):

    def __init__(self, num_nodes):
        sbeam.BeamModel.__init__(self)

        self.num_nodes = num_nodes
        self.num_elements = num_nodes - 1
        self.num_dofs = self.num_nodes * 3
        self.E = 2.1e11 #Youngs modulus N/m^2
        self.A = 0.1 #Area (unit width) m^2
        self.I = 1/12000 #Moment of inertia m^4
        self.ep = np.array([self.E, self.A, self.I])
        self.L = 1.6 # support distance m
        self.R = 1 # arch radius
        self.H = 0.4 # arch height

        theta_0 = math.asin((self.R-self.H)/self.R) #arch start angle, with respect to origin and horizontal
        theta_el = (math.pi - 2* theta_0)/self.num_elements #element angle increment



        self.coords = np.zeros((self.num_nodes,2)) # Coordinates for all the nodes
        self.dispState = np.zeros(self.num_nodes*3)
        self.Ndofs = np.zeros((self.num_nodes,3)) # Dofs for all the nodes

        for i in range(self.num_nodes):
            theta = theta_0 + i*theta_el
            self.coords[i,0] = self.L / 2 - self.R * math.cos(theta) #x-coordinates
            self.coords[i,1] = (self.H + self.R * (math.sin(theta) - 1) )  #y-coordinates
            self.Ndofs[i,:] = np.array([1,2,3],dtype=int) + i*3

        self.Edofs = np.zeros((self.num_elements,6),dtype=int) # Element dofs for all the elements
        self.Enods = np.zeros((self.num_elements,2),dtype=int) # Element nodes for all the elements
        for i in range(self.num_elements):
            self.Edofs[i,:] = np.array([1,2,3,4,5,6],dtype=int) + 3*i
            self.Enods[i,:] = np.array([1,2],dtype=int) + i

        # Fix x and y at first node and y at last node
        self.bc = np.array([1,2,(self.num_nodes*3 -1), (self.num_nodes*3 -2)],dtype=int)

        # The external incremental load (linear scaling with lambda)
        mid_node      = (self.num_nodes +1) // 2
        mid_y_dof_idx = (mid_node-1) * 3 + 1
        self.inc_load = np.zeros(self.num_dofs)
        self.inc_load[mid_y_dof_idx] = -500.0e+7
        self.plotDof = mid_y_dof_idx + 1

# ------------------ Perform linear solution

num_nodes = 5
#beamModel = sbeam.CantileverWithEndMoment(num_nodes)
beamModel = DeepArchModel(num_nodes)

load_steps=0.01 
N_steps=50 
max_iter=30

archLength=0.02
max_steps=50
max_iter=30
#sbeam.solveNonlinLoadControl(beamModel,load_steps, N_steps, max_iter)
sbeam.solveArchLength(beamModel, archLength, max_steps, max_iter)
#sbeam.solveLinearSteps(beamModel)

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
# ------------------ implementing non-linear solution