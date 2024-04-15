from ABD_Killian import Material, Lamina, Layup, Laminate, GenerateLaminate, ABD_Matrix, T
import numpy as np

np.set_printoptions(precision= 7)

def Stress(laminate, F):
    # calculate laminate strains
    e = np.linalg.solve(laminate.ABD, F) # midplane strains and curvature full
    'global strains'
    for i, lamina in enumerate(laminate.laminae):
        #global strains
        max_strain     = e[:3] + max([lamina.z0, lamina.z1], key=abs) * e[3:]
        min_strain     = e[:3] + min([lamina.z0, lamina.z1], key=abs) * e[3:]

        lamina.strains  = [max_strain, min_strain]
        # local stresses
        max_stress = T(lamina.orientation)[0] @ (lamina.Q @ lamina.strains[0])  # lamina s1, s2, s6
        min_stress = T(lamina.orientation)[0] @ (lamina.Q @ lamina.strains[1])

        lamina.stresses = [max_stress, min_stress] # max stress is used in failure envelope calculations 

        laminate.maxstrains[:, i][:, None]  = T(lamina.orientation)[1] @ lamina.strains[0]
        laminate.minstrains[:, i][:, None]  = T(lamina.orientation)[1] @ lamina.strains[1]  # lamina e1, e2, e6
        laminate.maxstresses[:,i][:, None]  = lamina.stresses[0]                             # storing stress states into laminate as well
        laminate.minstresses[:,i][:, None]  = lamina.stresses[1]                             # storing stress states into laminate as well
        
