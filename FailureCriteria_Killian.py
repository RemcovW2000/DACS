import numpy as np

np.set_printoptions(linewidth=300, precision = 3)

def FibreFailure(stress, material): # stress = np.array with s1,s2,s6 ; mat_properties with [Xt,Xc,Yt,Yc,S]
    print('Fibre Failure')
    Xt, Xc, Yt, Yc, S = material.strengths
    s1, s2, s6        = stress
    if s1 < 0:
        R = -Xc
    else:
        R = Xt
    # m_cfrp = 1.1, m_gfrp = 1.3
    f = 1/(R)*(s1 - (material.v21 - material.v21_fibre*1.1*material.E1/material.E1_fibre)*s2)
    return f

def MatrixFailure(stress, material):
    Yt, Yc, S = material.strengths[2:]
    s2, s6        = stress[1:]
    if s2 > 0: # mode A
        print('mode A')
        f = np.sqrt((s6/S)**2 + (1 - 0.3*Yt/S)**2*(s2/Yt)**2) + 0.3*s2/S
    
    s23A = S/(2*0.2)*(np.sqrt(1+2*0.2*Yc/S) - 1)
    p12_minus = 0.25
    p23_minus  = p12_minus*s23A/S # p12(-)????
    s12c = S*np.sqrt(1 + 2*p23_minus)

    if abs(s2/s6) <= s23A/abs(s12c) and abs(s2/s6) >= 0:
        print('mode B')
        f = (np.sqrt(s6**2 + (p12_minus*s2)**2) + p12_minus*s2)/S

    else:
        print('mode C')
        f = ((s6/(2*(1 + p23_minus*S)))**2 + (s2/Yc)**2)*(Yc/-s2)
    return f