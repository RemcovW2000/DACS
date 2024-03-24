import numpy as np
from Laminate import Laminate
from Lamina import Lamina

np.set_printoptions(linewidth=300, precision = 2)

#now we test the code:
E1 = 140000
E2 = 10000
G12 = 5000
v12 = 0.3

elasticproperties = [E1, E2, G12, v12]

E11f = 500000
v21f = 0.1
msf = 1.1
R11t = 1000
R11c = 800
yt = 100
yc = 100
S = 100

failureproperties = [E11f, v21f, msf, R11t, R11c, yt, yc, S]
s0 = Lamina(0.25, 0, elasticproperties, failureproperties)
s1 = Lamina(0.25, 45, elasticproperties, failureproperties)
s2 = Lamina(0.25, -45, elasticproperties, failureproperties)
s3 = Lamina(0.25, 0, elasticproperties, failureproperties)

# create the laminas list, for the laminate function:
laminas = [s0, s1, s2, s3]

# creating the laminate object:
laminate = Laminate(laminas)
# now we can apply loads to the laminate: (in Newton per meter or newtonmeter per meter)
loadingratio = np.array([[1], [0.1], [0], [0], [0], [0]])

print(laminate.ProgressiveDamageAnalysis(loadingratio, 0.1))
