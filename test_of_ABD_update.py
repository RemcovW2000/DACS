from ABD_Killian import Material, Lamina, Layup, Laminate, GenerateLaminate, ABD_Matrix
import numpy as np
#specify lay-up(including symmetry and thicknesses), materials
lay_up    = Layup([0,90], 1,  [0.125]) 
material  = Material('UD CFRP', 140e3, 10e3, 5e3, 0.3, None, [1.5e3, 1.2e3, 0.05e3, 0.25e3, 0.07e3], 0.2, 230e3) # MPa, [-]
print('lay_up:', lay_up.stack)

# initialize lamina and laminate
lamina   = Lamina(material, None, None, None, None, None, None, None, 0)
laminate = Laminate('cross-ply', lay_up.stack, None, np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), np.empty((6,6)), np.empty((6,6)), 0, 0, 0, 0, 0, np.empty((3,len(lay_up.stack))), np.empty((3,len(lay_up.stack))))

# assign lamina orientations and thicknesses based upon specified lay-up 
GenerateLaminate(lamina, lay_up, laminate) #, laminae)
# calculate ABD
ABD_Matrix(lay_up, laminate)

print(laminate.ABD)

lay_up    = Layup([90], 1,  [0.125]) 
material  = Material('UD CFRP', 140e3, 10e3, 5e3, 0.3, None, [1.5e3, 1.2e3, 0.05e3, 0.25e3, 0.07e3], 0.2, 230e3) # MPa, [-]
print('lay_up:', lay_up.stack)

# initialize lamina and laminate
lamina   = Lamina(material, None, None, None, None, None, None, None, 0)
laminate = Laminate('cross-ply', lay_up.stack, None, np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), np.empty((6,6)), np.empty((6,6)), 0, 0, 0, 0, 0, np.empty((3,len(lay_up.stack))), np.empty((3,len(lay_up.stack))))

# assign lamina orientations and thicknesses based upon specified lay-up 
GenerateLaminate(lamina, lay_up, laminate) #, laminae)
# calculate ABD
ABD_Matrix(lay_up, laminate)

print(laminate.ABD)