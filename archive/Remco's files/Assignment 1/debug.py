import numpy as np
from Laminate import Laminate
from Lamina import Lamina
import matplotlib.pyplot as plt

np.set_printoptions(linewidth=300, precision=2)

# now we test the code:
E1 = 145.3e3     # From assignment
E2 = 8.5e3       # From assignment
G12 = 4.58e3      # From assignment
v12 = 0.31      # From assignment

elasticproperties = [E1, E2, G12, v12]

# properties needed for failure analysis
E11f = 324e3   # TBD
v21f = 0.2      # TBD
msf = 1.1       # TBD
R11t = 1932     # From assignment
R11c = 1480     # From assignment
yt = 108        # From assignment
yc = 220        # From assignment
S = 132.8       # From assignment

failureproperties = [E11f, v21f, msf, R11t, R11c, yt, yc, S]
s0 = Lamina(0.25, 0, elasticproperties, failureproperties)
s1 = Lamina(0.25, 0, elasticproperties, failureproperties)
s2 = Lamina(0.25, 0, elasticproperties, failureproperties)
s3 = Lamina(0.25, 0, elasticproperties, failureproperties)

# create the laminas list, for the laminate function:
laminas = [s0, s1, s2, s3]

# creating the laminate object:
laminate = Laminate(laminas)
print(np.round(laminate.ABD_matrix, 2))

angle = 177
anglerad = np.deg2rad(angle)
loadingratio = np.array([[0],
                                     [np.sin(np.deg2rad(anglerad))],
                                     [np.cos(np.deg2rad(anglerad))],
                                     [0],
                                     [0],
                                     [0]])
laminate.progressive_damage_analysis(loadingratio, 0.1)

