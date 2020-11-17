# example Beam_models
# ----------------------------------------------------------------
# PURPOSE
#  Starting point for a couple of beam models


import math
import numpy as np
import matplotlib.pyplot as plt
import CorotBeam_with_TODO as CorotBeam
import matplotlib.animation as anm
from copy import deepcopy
# ----- Topology -------------------------------------------------

def solveArchLength(problem, archLength=0.02, max_steps=50, max_iter=30):
    num_dofs = problem.get_num_dofs()
    uVec = np.zeros(num_dofs)
    res_Vec = np.zeros(num_dofs)
    Lambda = 0.0

    d_q_prev = np.zeros(num_dofs)

    for iStep in range(max_steps):
        #TODO: Implement this predictor step, (sverre: trur eg er ferdig)
        q_Vec = problem.get_incremental_load(Lambda)
        K_mat = problem.get_K_sys(uVec)

        w_q0 = np.linalg.solve(K_mat,q_Vec)
        f = math.sqrt(1 + w_q0.T @ w_q0)

        if (w_q0.T @ uVec) > 1:
            delta_Lambda = archLength / f
        else:
            delta_Lambda = - archLength / f

        Lambda += delta_Lambda
        uVec += (delta_Lambda * w_q0)

        for iIter in range(max_iter):
            res_Vec = problem.get_residual(Lambda, uVec)
            K_mat = problem.get_K_sys(uVec)
            q_Vec = problem.get_incremental_load(Lambda)

            w_q = np.linalg.solve(K_mat,q_Vec)
            w_r = np.linalg.solve(K_mat,-res_Vec)


            d_Lambda = - w_q.T @ w_r / (1 + w_q.T @ w_r)

            Lambda += d_Lambda
            uVec += w_r + d_Lambda*w_q

            # TODO: Implement this corrector step, (sverre: trur eg er ferdig, men funker ikkje endo)

            res_Vec = problem.get_residual(Lambda , uVec)
            if (res_Vec.dot(res_Vec) < 1.0e-15):  #check if residual is small enough
                break

        problem.append_solution(Lambda, uVec)
        print(" ")

def solveNonlinLoadControl(problem, load_steps=0.01, max_steps=100, max_iter=30):
    num_dofs = problem.get_num_dofs()
    uVec   = np.zeros(shape=(num_dofs,))
    #d_uVec = np.zeros(shape=(num_dofs,1))

    for iStep in range(max_steps):
        
        #TODO: Implement this Predictor: (Almar: Load control, forward euler)

        Lambda = load_steps * iStep
        Lambda_nxt = load_steps * (iStep+1)
        q_Vec   = problem.get_incremental_load(Lambda)
        K_mat = problem.get_K_sys(uVec)

        Delta_Lambda = Lambda_nxt-Lambda
        K_mat_inv = np.linalg.inv(K_mat)
        v_mat = K_mat_inv @ q_Vec
        d_uVec = v_mat * Delta_Lambda

        uVec = uVec + d_uVec
        
        for iIter in range(max_iter):

            # TODO: Implement this Korrektor (Almar: Newton metode her)
            # Husk at load  control betyr: at all kerreksjon er horisontal (lasten holdes konstant i koreksjonen)
            # Her regnes residual derivasjonen altså ikke eksakt, men tilnermes med en delta i load_steps, sikkert ikke like bra
            

            res_Vec = problem.get_residual(Lambda, uVec)
            if (res_Vec.dot(res_Vec) < 1.0e-15):
                break 
            #d_res_Vec = (problem.get_residual(Lambda,uVec+load_steps) - res_Vec)/load_steps #Finner endring i residual. (dette er litt fishi)
            delta_uVec = K_mat_inv @ (-res_Vec)
            uVec = uVec + delta_uVec
            
            K_mat = problem.get_K_sys(uVec)
            K_mat_inv = np.linalg.inv(K_mat) 
            

        problem.append_solution(Lambda, uVec)
        print(" ")

        problem.append_solution(Lambda, uVec)
        print("Non-Linear load step {:}  load_factor= {:12.3e}".format(iStep, Lambda))



def solveLinearSteps(problem, load_steps=0.01, max_steps=100):
    num_dofs = problem.get_num_dofs()
    uVec = np.zeros(num_dofs)

    for iStep in range(max_steps):

        Lambda = load_steps * iStep

        q_Vec   = problem.get_incremental_load(Lambda)

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
        K_sys = np.zeros((self.num_dofs,self.num_dofs))

        for iel in range(self.num_elements):
            inod1 = self.Enods[iel,0]-1
            inod2 = self.Enods[iel,1]-1
            ex1 = self.coords[inod1,0]
            ex2 = self.coords[inod2,0]
            ex = np.array([ex1,ex2])
            ey = np.array([self.coords[inod1,1],self.coords[inod2,1]])
            Edofs = self.Edofs[iel] - 1
            Ke = CorotBeam.beam2e(ex, ey, self.ep)
            #Ke = CorotBeam.beam2corot_Ke_and_Fe(ex, ey, self.ep, disp_sys[np.ix_(Edofs)] )[0] #korotert stivhet
            K_sys[np.ix_(Edofs,Edofs)] += Ke

        # Set boundary conditions
        for idof in range(len(self.bc)):
            idx = self.bc[idof] - 1
            K_sys[idx,:]   = 0.0
            K_sys[:,idx]   = 0.0
            K_sys[idx,idx] = 1.0

        return K_sys

    def get_K_sys(self, disp_sys):
        # Build system stiffness matrix for the structure
        K_sys = np.zeros((self.num_dofs,self.num_dofs))

        for iel in range(self.num_elements):
            inod1 = self.Enods[iel,0]-1 
            inod2 = self.Enods[iel,1]-1
            ex1 = self.coords[inod1,0]
            ex2 = self.coords[inod2,0]
            ex = np.array([ex1,ex2])
            ey = np.array([self.coords[inod1,1],self.coords[inod2,1]])
            Edofs = self.Edofs[iel] - 1
            #Ke = CorotBeam.beam2e(ex, ey, self.ep)
            Ke = CorotBeam.beam2corot_Ke_and_Fe(ex, ey, self.ep, disp_sys[np.ix_(Edofs)] )[0] #korotert stivhet
            K_sys[np.ix_(Edofs,Edofs)] += Ke

        # Set boundary conditions
        for idof in range(len(self.bc)):
            idx = self.bc[idof] - 1
            K_sys[idx,:]   = 0.0
            K_sys[:,idx]   = 0.0
            K_sys[idx,idx] = 1.0

        return K_sys

    def get_num_dofs(self):
        num_dofs = self.num_nodes * 3
        return num_dofs

    def get_internal_forces(self, disp_sys):
        # Build internal force vector for the structure
        f_int_sys = np.zeros(self.num_dofs)

        for iel in range(self.num_elements):
            inod1 = self.Enods[iel,0]-1
            inod2 = self.Enods[iel,1]-1
            ex1 = self.coords[inod1,0]
            ex2 = self.coords[inod2,0]
            ex = np.array([ex1,ex2])
            ey = np.array([self.coords[inod1,1],self.coords[inod2,1]])

            Edofs = self.Edofs[iel] - 1 #Tror ikke vi har forstaatt Edofs helt. i=1 -> Edofs = [0,1,2,3,4,5], i=2 -> Edofs = [3,4,5,6,7,8]
            f_int_e = CorotBeam.beam2corot_Ke_and_Fe(ex, ey, self.ep, disp_sys[np.ix_(Edofs)])[1] #korotert internal force [6x1]
            #disp_e = disp_sys[np.ix_(Edofs)] 
            #f_int_e = Ke * disp_e
            f_int_sys[np.ix_(Edofs)] += f_int_e #Kan hende dette maa endres [num_dofs,1]

        return f_int_sys

    def get_incremental_load(self,loadFactor):
        return self.inc_load

    def get_external_load(self,loadFactor):
        return (self.inc_load * loadFactor)

    def get_residual(self,loadFactor,disp_sys):
        f_int = self.get_internal_forces(disp_sys)
        f_res = self.get_external_load(loadFactor) + f_int
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

        ax.plot(plottedDeformations,self.load_history)
        ax_shape.plot(x,y)

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
        #self.inc_load[-1] = 1.0e6
        self.inc_load[-1] = MomentFullCircle
        self.plotDof = self.num_dofs - 1 # Setting which dof for the Load-disp curve: y disp of last node


'''
# ------------------ Perform linear solution

num_nodes = 9
beamModel = SimplySupportedBeamModel(num_nodes)
#beamModel = CantileverWithEndMoment(num_nodes)

solveLinearSteps(beamModel)

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
