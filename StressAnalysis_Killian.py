from ABD_Killian import T
import numpy as np
from uncertainties import ufloat

def Stress(laminate, F):
    # calculate laminate strains
    e = np.linalg.solve(laminate.ABD, F) # midplane strains and curvature full 
    for i, lamina in enumerate(laminate.laminae):
        lamina.strains  = e[:3] + max([lamina.z0, lamina.z1], key=abs) * e[3:]    # lamina e1, e2, e6
        lamina.stresses = T(lamina.orientation)[0] @ (lamina.Q @ lamina.strains)  # lamina s1, s2, s6
        #print(lamina.Q)
        laminate.strains[:, i][:, None]  = lamina.strains
        laminate.stresses[:,i][:, None]  = lamina.stresses                         # storing stress states into laminate as well