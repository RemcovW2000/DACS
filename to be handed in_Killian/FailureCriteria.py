import numpy as np
from uncertainties import ufloat

np.set_printoptions(linewidth=300, precision = 3)

def MaxStressFailure(stress, material): # stress = np.array with s1,s2,s6 ; mat_properties with [Xt,Xc,Yt,Yc,S]
    Xt, Xc, Yt, Yc, S = material.strengths
    s1, s2, s6        = stress

    if s1 < 0:
        FI1 = abs(s1)/Xc # failure in X corresponds to fibre failure
    else:
        FI1 = abs(s1)/Xt

    if s2 < 0:
        FI2 = abs(s2)/Yc # failure in Y corresponds to matrix failure
    else:
        FI2 = abs(s2)/Yt

    FI3 = abs(s6)/abs(S)

    FI = {'FI1':FI1[0], 'FI2':FI2[0], 'FI3':FI3[0]} 

    return FI

def PuckFailure(stress, material, sanity_check):
    FI1, FI2, FI3, FI4 = np.array([0]), np.array([0]), np.array([0]), np.array([0])
    Xt, Xc, Yt, Yc, S = material.strengths
    s1, s2, s6        = stress

    p12_minus  = 0.2
    s23A       = S/(2*p12_minus)*(np.sqrt(1 + 2*p12_minus*Yc/S) - 1)
    p23_minus  = p12_minus*s23A/S
    s12c       = S*np.sqrt(1 + 2*p23_minus)
    #theta_fp   = np.rad2deg(np.arccos(np.sqrt(s23A/-s2)))

    'FIBRE FAILURE'
    if s1 < 0:
        R = -Xc
    else:
        R = Xt
    # m_cfrp = 1.1, m_gfrp = 1.3 #TODO: check material assignment
    FI1 = 1/(R)*(s1 - (material.v21 - material.v21_fibre*1.3*material.E1/material.E1_fibre)*s2)
    'MATRIX FAILURES: A, B or C'
    if s2 >= 0: # mode A
        if sanity_check:
            print('mode A')
        FI2 = np.sqrt((s6/S)**2 + (1 - 0.3*Yt/S)**2*(s2/Yt)**2) + 0.3*s2/S
    
    elif abs(s2 / (abs(s6) + 1e-11)) <= s23A / abs(s12c) and abs(s2 / (abs(s6) + 1e-11)) >= 0:
        if sanity_check:
            print('mode B')
        FI3 = (np.sqrt(s6**2 + (p12_minus*s2)**2) + p12_minus*s2) / S

    else:
        if sanity_check:
            print('mode C')
        FI4 = (((s6/(2*(1 + p23_minus)*S)))**2 + (s2/Yc)**2) * (Yc/-s2)

    FI = {'FI1':FI1[0], 'FI2':FI2[0],  'FI3':FI3[0],'FI4': FI4[0]} 

    return FI


def DegradationRule(lamina):
    value = 1e-20
    material = lamina.material
    if lamina.failure_counter != 0:
        if lamina.failure_counter == 1:
            if lamina.failure_mode == 'FI1': # same label holds for Puck and Max Stress Fiber Failure
                material.E1  *= value
                material.E2  *= value
                material.G12 *= value
                material.v12 *= value
                material.v21 *= value
            else:
                material.E2  *= 0.05
                material.v21 *= 0.05
                material.G12 *= 0.05 # most damages start in matrix, shear is determined by matrix, hence modulus also degraded
    
        if lamina.failure_counter == 2:
            material.E1  *= value
            material.E2  *= value
            material.G12 *= value
            material.v12 *= value
            material.v21 *= value





