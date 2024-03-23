from ABD_Killian import T
import numpy as np

def Stress(laminate, F):
    # calculate laminate strains
    e = np.linalg.solve(laminate.ABD, F) # midplane strains and curvature full 
    for lamina in laminate.laminae:
        lamina.strains = e[:3] + max([lamina.z0, lamina.z1], key=abs) * e[3:]    # lamina e1, e2, e6
        lamina.stresses = T(lamina.orientation)[0] @ (lamina.Q @ lamina.strains) # lamina s1, s2, s6
