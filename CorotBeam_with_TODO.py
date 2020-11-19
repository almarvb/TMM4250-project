# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 16:43:51 2018

@author: bjohau
"""

import numpy as np
import math

def rot_matrix(theta):
    """
    Return the 2x2 rotation matrix representing a rotation theta
    :param theta:  rotation angle in radians
    :return: Rotation matrix (or tensor)
    """
    s = math.sin(theta)
    c = math.cos(theta)
    R = np.array([[c, -s],
                  [s,  c]])
    return R

def beam2local_def_disp(ex,ey, disp_global):
    """

    :param ex: element x coordinate [x1, x2] in undeformed position
    :param ey: element y coordinate [y1, y2] in undeformed position
    :param disp_global:  displacement vector [u1, v1, r1, u2, v2, r2] in global directions
    :return: disp_local_def: displacement vector [u1, v1, r1, u2, v2, r2] in local directions
    """
    eVec12 = np.array([ex[1] - ex[0], ey[1] - ey[0]])
    L0 = math.sqrt(eVec12 @ eVec12)

    # TODO: Quite a bit here (Almar: tror jeg har gjort alt)
    ex0 = np.array([[(ex[1]-ex[0])/L0],
                    [(ey[1]-ey[0])/L0]])
    exn = np.array([[(ex[1]+disp_global[3])-(ex[0]+disp_global[0])],
                    [(ey[1]+disp_global[4])-(ey[0]+disp_global[1])]])

    Ld = math.sqrt(exn.T @ exn)
    exn /= Ld

    eyn = np.array([[-exn[1,0]],
                    [exn[0,0]]])
    R1 = rot_matrix(disp_global[2])
    R2 = rot_matrix(disp_global[5])
    t1 = R1 @ ex0
    t2 = R2 @ ex0
    theta1_def = math.asin(eyn.T @ t1)
    theta2_def = math.asin(eyn.T @ t2)

    def_disp_local = np.array([ -0.5*(Ld - L0),
                                0.0,
                                theta1_def,
                                0.5 * (Ld - L0),
                                0.0,
                                theta2_def])
    return def_disp_local


def beam2corot_Ke_and_Fe(ex,ey,ep, disp_global):
    """
    Compute the stiffness matrix and internal forces for a two dimensional beam element
    relative to deformed configuration.
    
    :param list ex: element x coordinates [x1, x2]
    :param list ey: element y coordinates [y1, y2]
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia
    :param list disp_global displacement vector for the element [ux1,uy1,rz1,ux2,uy2,rz2]


    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: element internal force vector [6 x 1]
    """
    # TODO: Quite a bit here (Sverre: prøver å fnne ut av denne. vanskelig)
    # Undeformed length and unit vector along element
    eVec12 = np.array([ex[1] - ex[0], ey[1] - ey[0]])
    L0 = math.sqrt(eVec12 @ eVec12)
    eVec12 /= L0 #Element undeformed length

    # Deformed position, unit vector and length along element
    ex_def = ex + [disp_global[0], disp_global[3]]
    ey_def = ey + [disp_global[1], disp_global[4]]

    eVec12_def = np.array([ex_def[1] - ex_def[0], ey_def[1] - ey_def[0]])
    Ld = math.sqrt(eVec12_def @ eVec12_def) #Deformed element length



    disp_def_local = beam2local_def_disp(ex,ey,disp_global)

    Kle = beam2local_stiff(L0,ep) # Element material stiffness of undeformed ghost element in local



    f_int_lin = Kle @ disp_def_local #find internal forces from linear stiffness
    N = f_int_lin[3] #pick out normal force and shear force to build geometric stiffness matrix
    V = f_int_lin[4]  #Nå oppdateres kun V for hver iterasjon, N blir konstant og kjempe stor!!!! <--- SE På DETTE SVERRE

    Kg_sym = np.array([
                        [ 0        , -V/(2*Ld), 0 , 0       , V/(2*Ld) , 0 ],
                        [ -V/(2*Ld), N/Ld     , 0 , V/(2*Ld), -N/Ld    , 0 ],
                        [ 0        , 0        , 0 , 0       , 0        , 0 ],
                        [ 0        , V/(2*Ld) , 0 , 0       , -V/(2*Ld), 0 ],
                        [ V/(2*Ld) , -N/Ld    , 0 ,-V/(2*Ld), N/Ld     , 0 ],
                        [0         , 0        , 0 , 0       , 0        , 0 ]
                        ]) #build element symmetric geometric stiffness matrix

    Te = beam2corot_Te(ex_def,ey_def) #Transformation matrix

    #Ke_g = Te.T @ Kg_sym @ Te #geometric stiffness, global coordinates
    #Ke_m = Te.T @ Kle @ Te #material stiffness, global coordinates

    #Ke_global =  Ke_m  + Ke_g   #element stiffness, global coordinates
    K_loc = Kle + Kg_sym
    #K_loc = Kle   #TODO, delete this !!!
    Ke_global =  Te.T @ K_loc @ Te
    #Ke_global = Te.T @ Kle @ Te  # Tar kun med matrial stivhet
    fe_int_global = Te.T @ f_int_lin #Internal forces, global coordinates

    return Ke_global, fe_int_global

    
def beam2corot_Te(ex,ey):
    """
    Compute the transformation matrix for an element
    
    :param list ex: element x coordinates [x1, x2]
    :param list ey: element y coordinates [y1, y2]
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia   
    :param list eq: distributed loads, local directions [qx, qy]
    :return mat Te: element transformation from global to local
    """

    n = np.array([ex[1]-ex[0],ey[1]-ey[0]])
    L = np.linalg.norm(n)
    n = n / L  
    
    Te=np.array([
        [ n[0], n[1],  0.,    0.,    0.,   0.],
        [-n[1], n[0],  0.,    0.,    0.,   0.],
        [0.,    0.,    1.,    0.,    0.,   0.],
        [0.,    0.,    0.,   n[0],  n[1],  0.],
        [0.,    0.,    0.,  -n[1],  n[0],  0.],
        [0.,    0.,    0.,    0.,    0.,   1.]
    ])
    

    return Te
    
    
def beam2local_stiff(L,ep):
    """
    Compute the stiffness matrix for a two dimensional beam element.
    
    :param list L : element length
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia   
    :return mat Kle: element stiffness matrix [6 x 6]
    """
    
    E=ep[0]
    A=ep[1]
    I=ep[2]
        
    Kle = np.array([
        [E*A/L,              0.,           0.,    -E*A/L,            0.,           0. ],
        [   0.,    12*E*I/L**3.,  6*E*I/L**2.,        0., -12*E*I/L**3.,  6*E*I/L**2. ],
        [   0.,     6*E*I/L**2.,      4*E*I/L,        0.,  -6*E*I/L**2.,     2*E*I/L  ],
        [-E*A/L,             0.,           0.,     E*A/L,            0.,           0. ],
        [   0.,   -12*E*I/L**3., -6*E*I/L**2.,        0.,  12*E*I/L**3., -6*E*I/L**2. ],
        [   0.,     6*E*I/L**2.,      2*E*I/L,        0.,  -6*E*I/L**2.,      4*E*I/L ]
    ])
     
    return Kle


def beam2e(ex, ey, ep, eq=None):
    """
    Compute the linear stiffness matrix for a two dimensional beam element.
    Largely from CALFEM core module

    :param list ex: element x coordinates [x1, x2]
    :param list ey: element y coordinates [y1, y2]
    :param list ep: element properties [E, A, I], E - Young's modulus, A - Cross section area, I - Moment of inertia
    :param list eq: distributed loads, local directions [qx, qy]
    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: element stiffness matrix [6 x 1] (if eq!=None)
    """

    n = np.array([ex[1] - ex[0], ey[1] - ey[0]])
    L = np.linalg.norm(n)
    n = n / L

    qx = 0.
    qy = 0.
    if not eq is None:
        qx = eq[0]
        qy = eq[1]

    Kle = beam2local_stiff(L,ep)

    fle = L * np.mat([qx / 2, qy / 2, qy * L / 12, qx / 2, qy / 2, -qy * L / 12]).T

    Te = beam2corot_Te(ex,ey)

    Ke = Te.T @ Kle @ Te
    fe = Te.T @ fle

    if eq is None:
        return Ke
    else:
        return Ke, fe