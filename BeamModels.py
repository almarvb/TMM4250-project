'''This file contains the BeamModel class with three beam models:
 CantileverWithEndMoment, SimplySupportedBeamModel and DeepArchModel,
 along with the three solving algotithms:
 SolveLinearSteps, SolveNonlinLoadControl and SolveArcLength.'''


import math
import numpy as np
import matplotlib.pyplot as plt
import CorotBeam as CorotBeam
import matplotlib.animation as anm
from copy import deepcopy
# ----- Topology -------------------------------------------------


def solveArcLength(problem, archLength, max_steps, max_iter):
    num_dofs = problem.get_num_dofs()
    uVec = np.zeros(num_dofs)
    res_Vec = np.zeros(num_dofs)
    Lambda = 0.0

    v_0 = np.zeros(num_dofs)
    q_Vec = problem.get_incremental_load(Lambda)
    for iStep in range(max_steps):  # Predictor step

        K_mat = problem.get_K_sys(uVec)

        w_q0 = np.linalg.solve(K_mat, q_Vec)
        f = math.sqrt(1 + w_q0.T @ w_q0)

        if (w_q0.T @ v_0) >= 0.0:
            delta_Lambda = archLength / f
        else:
            delta_Lambda = - archLength / f

        v_0 = (delta_Lambda * w_q0)
        Lambda += delta_Lambda
        uVec += v_0

        bConverged = False
        res_Vec = problem.get_residual(Lambda, uVec)

        for iIter in range(max_iter):  # Corrector step
            K_mat = problem.get_K_sys(uVec)  # Build system stiffness matrix

            w_q = np.linalg.solve(K_mat, q_Vec)
            w_r = np.linalg.solve(K_mat, -res_Vec)

            d_Lambda = - (w_q.T @ w_r) / (1 + w_q.T @ w_q)

            Lambda += d_Lambda
            d_uVec = w_r + d_Lambda*w_q
            uVec += (d_uVec)

            res_Vec = problem.get_residual(Lambda, uVec)
            resNorm = res_Vec.dot(res_Vec)
            print("iter {:}  resNorm= {:12.3e}".format(iIter, resNorm))
            if (resNorm < 1.0e-15):
                bConverged = True  # check if residual is small enough
                break
        if (not bConverged):
            print("Did not converge !!!!!!!!")
            break
        else:
            problem.append_solution(Lambda, deepcopy(uVec))
            print("Arc Length step {:}  load_factor= {:12.3e}".format(iStep, Lambda))

def solveNonlinLoadControl(problem, load_steps, max_steps, max_iter):
    num_dofs = problem.get_num_dofs()
    
    uVec   = np.zeros(num_dofs)
    
    for iStep in range(max_steps):
        
        #Increases the load factor Lambda for each iStep

        Lambda = load_steps * iStep

        bConverged = False
        
        for iIter in range(max_iter):
            #Keeps Lambda constant while iterting ftowards equilibrium
            res_Vec = problem.get_residual(Lambda, uVec)

            resNorm = res_Vec.dot(res_Vec)
            print("iter {:}  resNorm= {:12.3e}".format(iIter, resNorm))  
            if (resNorm < 1.0e-8):
                bConverged = True
                break 
            
            K_mat = problem.get_K_sys(uVec)

            delta_uVec = np.linalg.solve(K_mat, res_Vec)
            uVec = uVec + delta_uVec


        if (not bConverged):
            print("Did not converge !!!!!!!!")
            break
        else:
            problem.append_solution(Lambda, deepcopy(uVec))
            print("Non-Linear load step {:}  load_factor= {:12.3e}".format(iStep, Lambda))

def solveLinearSteps(problem, load_steps, max_steps):
    num_dofs = problem.get_num_dofs()
    uVec = np.zeros(num_dofs)

    for iStep in range(max_steps):

        Lambda = load_steps * iStep

        q_Vec = problem.get_incremental_load(Lambda)

        K_mat = problem.get_K_sys_lin(uVec)

        d_q = np.linalg.solve(K_mat, q_Vec)

        uVec = d_q * Lambda

        problem.append_solution(Lambda, uVec)
        print("Linear load step {:}  load_factor= {:12.3e}".format(iStep, Lambda))


class BeamModel:

    def __init__(self):
        self.inc_load = None
        self.bc = None        # List of fixed dofs, 1 based
        self.num_nodes = None
        self.coords = None    # Nodal coordinates, array (num_nodes x 2)
        self.Edofs = None     # Element dofs, 1 based    (num_elements x 6)
        self.Enods = None     # Element nodes, 1 based   (num_elements x 2)
        self.ep = None        # Element properties, array of 3 values
        self.num_elements = None
        self.num_dofs = None

        # Plotting related data
        self.plotDof = None   # Plotting dof for load displacment plot, 1 based
        self.plotScaleFactor = 1.0
        self.load_history = []
        self.disp_history = []

    def get_K_sys_lin(self, disp_sys):
        # Build system stiffness matrix for the structure
        K_sys = np.zeros((self.num_dofs, self.num_dofs))

        for iel in range(self.num_elements):
            inod1 = self.Enods[iel, 0]-1
            inod2 = self.Enods[iel, 1]-1
            ex1 = self.coords[inod1, 0]
            ex2 = self.coords[inod2, 0]
            ex = np.array([ex1, ex2])
            ey = np.array([self.coords[inod1, 1], self.coords[inod2, 1]])
            Edofs = self.Edofs[iel] - 1
            
            Ke = CorotBeam.beam2e(ex, ey, self.ep)
            K_sys[np.ix_(Edofs,Edofs)] += Ke

        # Set boundary conditions
        for idof in range(len(self.bc)):
            idx = self.bc[idof] - 1
            K_sys[idx, :]   = 0.0
            K_sys[:, idx]   = 0.0
            K_sys[idx, idx] = 1.0

        return K_sys

    def get_K_sys(self, disp_sys):
        # Build system stiffness matrix for the structure
        K_sys = np.zeros((self.num_dofs, self.num_dofs))

        for iel in range(self.num_elements):
            inod1 = self.Enods[iel, 0]-1
            inod2 = self.Enods[iel, 1]-1
            ex1 = self.coords[inod1, 0]
            ex2 = self.coords[inod2, 0]
            ex = np.array([ex1, ex2])
            ey = np.array([self.coords[inod1, 1], self.coords[inod2, 1]])
            Edofs = self.Edofs[iel] - 1

            Ke = CorotBeam.beam2corot_Ke_and_Fe(ex, ey, self.ep, disp_sys[np.ix_(Edofs)])[0]  # Stiffness of corotated element
            K_sys[np.ix_(Edofs, Edofs)] += Ke

        # Set boundary conditions
        for idof in range(len(self.bc)):
            idx = self.bc[idof] - 1
            K_sys[idx, :]   = 0.0
            K_sys[:, idx]   = 0.0
            K_sys[idx, idx] = 1.0

        return K_sys

    def get_num_dofs(self):
        num_dofs = self.num_nodes * 3
        return num_dofs

    def get_internal_forces(self, disp_sys):
        # Build internal force vector for the structure
        f_int_sys = np.zeros(self.num_dofs)

        for iel in range(self.num_elements):
            inod1 = self.Enods[iel, 0]-1
            inod2 = self.Enods[iel, 1]-1
            ex1 = self.coords[inod1, 0]
            ex2 = self.coords[inod2, 0]
            ex = np.array([ex1, ex2])
            ey = np.array([self.coords[inod1, 1], self.coords[inod2, 1]])

            Edofs = self.Edofs[iel] - 1 
            f_int_e = CorotBeam.beam2corot_Ke_and_Fe(ex, ey, self.ep, disp_sys[np.ix_(Edofs)])[1] #korotert internal force [6x1]

            f_int_sys[np.ix_(Edofs)] += f_int_e 

        return f_int_sys

    def get_incremental_load(self, loadFactor):
        return self.inc_load

    def get_external_load(self, loadFactor):
        return (self.inc_load * loadFactor)

    def get_residual(self, loadFactor, disp_sys):
        f_int = self.get_internal_forces(disp_sys)
        f_ext = self.get_external_load(loadFactor)
        f_res = f_ext - f_int
        
        # Set boundary conditions
        for idof in range(len(self.bc)):
            idx = self.bc[idof] - 1
            f_res[idx]   = 0.0 #For any fixed DOFs the residual must be set to 0

        return f_res

    def append_solution(self, loadFactor, disp_sys):
        self.load_history.append(loadFactor)
        self.disp_history.append(disp_sys)

    def plotDispState(self, step, limits=None, scaleFactor=1 ):
        # DOF: DOF number of displacement to be plotted
        # Scale: scale factor for displacement
        plt.close('all')  # Close all currently open figures
        # Getting deformations from all nodes:
        # deformations = self.disp_history[:,DOF]
        num_steps = len(self.load_history)
        num_nodes = self.num_nodes

        # Starting subplot:
        fig, (ax, ax_shape) = plt.subplots(nrows=1, ncols=2, figsize=(20, 10))

        # Setting up force-displacement curve:
        plottedDeformations = np.zeros(num_steps)
        for i in range(num_steps):
            plottedDeformations[i] = self.disp_history[i][self.plotDof-1]
        plottedDeformations -= plottedDeformations[0]
        plottedDeformations *= self.plotScaleFactor

        ax.plot(plottedDeformations, self.load_history, label="Node")
        ax.plot([plottedDeformations[step]], [self.load_history[step]], '.', color='red', markerSize=20)
        ax.legend(loc='best', fontsize=12)
        ax.set_xlabel('Displacement', fontsize=12)
        ax.set_ylabel('Applied force', fontsize=12)
        ax.set_title('Force-displacement')

        # Animating:
        x = np.zeros(num_nodes)
        y = np.zeros(num_nodes)
        if step >= num_steps:
            step = num_steps -1

        for i in range(num_nodes):
            x[i] = self.disp_history[step][0 + i * 3] + self.coords[i, 0]
            y[i] = self.disp_history[step][1 + i * 3] + self.coords[i, 1]
        line, = ax_shape.plot(x, y)
        if limits is None:
            ax_shape.axis('equal')
        else:
            ax_shape.set_xlim(limits[0], limits[1])
            ax_shape.set_ylim(limits[2], limits[3])

        ax.plot(plottedDeformations, self.load_history)
        ax_shape.plot(x, y)

        plt.show(block=True)

class SimplySupportedBeamModel(BeamModel):

    def __init__(self, num_nodes):
        BeamModel.__init__(self)

        self.num_nodes = num_nodes
        self.num_elements = num_nodes - 1
        self.num_dofs = self.num_nodes * 3
        self.E = 2.1e11
        self.A = 45.3e-4
        self.I = 2510e-8
        self.ep = np.array([self.E, self.A, self.I])
        self.L_total = 9.0

        L_el = self.L_total / self.num_elements

        self.coords = np.zeros((self.num_nodes,2)) # Coordinates for all the nodes
        self.dispState = np.zeros(self.num_nodes*3)
        self.Ndofs  = np.zeros((self.num_nodes,3)) # Dofs for all the nodes

        for i in range(self.num_nodes):
            self.coords[i,0] = i * L_el
            self.Ndofs[i,:] = np.array([1,2,3],dtype=int) + i*3

        self.Edofs = np.zeros((self.num_elements,6),dtype=int) # Element dofs for all the elements
        self.Enods = np.zeros((self.num_elements,2),dtype=int) # Element nodes for all the elements
        for i in range(self.num_elements):
            self.Edofs[i,:] = np.array([1,2,3,4,5,6],dtype=int) + 3*i
            self.Enods[i,:] = np.array([1,2],dtype=int) + i

        # Fix x and y at first node and y at last node
        self.bc = np.array([1,2,(self.num_nodes*3 -1)],dtype=int)

        # The external incremental load (linear scaling with lambda)
        mid_node      = (self.num_nodes +1) // 2
        mid_y_dof_idx = (mid_node-1) * 3 + 1
        self.inc_load = np.zeros(self.num_dofs)
        self.inc_load[mid_y_dof_idx] = 500.0e+4
        self.plotDof = mid_y_dof_idx + 1

class CantileverWithEndMoment(BeamModel):

    def __init__(self, num_nodes):
        BeamModel.__init__(self)

        self.num_nodes = num_nodes
        self.num_elements = num_nodes - 1
        self.num_dofs = self.num_nodes * 3
        self.E = 2.1e11
        self.A = 45.3e-4
        self.I = 2510e-8
        self.ep = np.array([self.E, self.A, self.I])
        self.L_total = 9.0

        MomentFullCircle = 2.0 * math.pi * self.E * self.I / self.L_total

        L_el = self.L_total / self.num_elements

        self.coords = np.zeros((self.num_nodes,2))   # Coordinates for all the nodes
        self.dispState = np.zeros(self.num_nodes*3)  # Current displacement state for all dofs
        self.Ndofs  = np.zeros((self.num_nodes,3))   # Dofs for all the nodes

        for i in range(self.num_nodes):
            self.coords[i,0] = i * L_el              # Setting the x-coordinat value
            self.Ndofs[i,:] = np.array([1,2,3],dtype=int) + i*3

        self.Edofs = np.zeros((self.num_elements,6),dtype=int) # Element dofs for all the elements
        self.Enods = np.zeros((self.num_elements,2),dtype=int) # Element nodes for all the elements
        for i in range(self.num_elements):
            self.Edofs[i,:] = np.array([1,2,3,4,5,6],dtype=int) + 3*i
            self.Enods[i,:] = np.array([1,2],dtype=int) + i

        # Fix x, y and rotation at first node
        self.bc = np.array([1,2,3],dtype=int)

        # The external incremental load (linear scaling with lambda)
        self.inc_load = np.zeros(self.num_dofs)
        #self.inc_load[-1] = 1.0e8
        self.inc_load[-1] = MomentFullCircle
        self.plotDof = self.num_dofs - 1 # Setting which dof for the Load-disp curve: y disp of last node

class DeepArchModel(BeamModel):

    def __init__(self, num_nodes):
        BeamModel.__init__(self)

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
