from ABD_Killian import Material, Lamina, Layup, Laminate, GenerateLaminate, ABD_Matrix #, T
from StressAnalysis_Killian import Stress
from FailureCriteria_Killian import FibreFailure, MatrixFailure
import numpy as np
import copy

## example analysis 

# specify lay-up(including symmetry and thicknesses), materials
lay_up    = Layup([0,45,-45,90], 1,  [0.2]) 
material  = Material('UD CFRP', 140e3, 10e3, 5e3, 0.3, None, [1.5e3, 1.2e3, 0.05e3, 0.25e3, 0.07e3], 0.2, 230e3) # GPa, [-]
print('lay_up:', lay_up.stack)

# initialize lamina and laminate
lamina   = Lamina(material, None, None, None, None, None, None, None)
laminate = Laminate('cross-ply', lay_up.stack, None, np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), np.empty((6,6)), np.empty((6,6)), 0, 0, 0, 0, 0)

# assign lamina orientations and thicknesses based upon specified lay-up 
GenerateLaminate(lamina, lay_up, laminate) #, laminae)
# calculate ABD
ABD_Matrix(lay_up, laminate)

# initialize force vector F, containing force and moment intensities
F = np.ones((6,1)) # (Nx, Ny, Ns, Mx, My, Ms)
# apply loadcase
F[0] = 10 # N/mm -> results in stresses in MPa (N/mm^2)

Stress(laminate, F)

for lamina in laminate.laminae:
    #print(lamina.stresses)
    print('##############################3')
    print(FibreFailure(lamina.stresses, material))
    print(MatrixFailure(lamina.stresses, material))