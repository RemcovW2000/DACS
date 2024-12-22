import numpy as np
from Laminate import Laminate
from Lamina import Lamina
from stringers import TStringer

np.set_printoptions(linewidth=300, precision=2)

# now we test the code:
E1 = 145300     # From assignment
E2 = 8500       # From assignment
G12 = 4580      # From assignment
v12 = 0.31      # From assignment

elasticproperties = [E1, E2, G12, v12]

# properties needed for failure analysis
E11f = 500000   # TBD
v21f = 0.1      # TBD
msf = 1.1       # TBD
R11t = 1932     # From assignment
R11c = 1480     # From assignment
yt = 108        # From assignment
yc = 220        # From assignment
S = 132.8       # From assignment

failureproperties = [E11f, v21f, msf, R11t, R11c, yt, yc, S]
V0 = Lamina(0.25, 0, elasticproperties, failureproperties)
V1 = Lamina(0.25, 0, elasticproperties, failureproperties)
V2 = Lamina(0.25, 0, elasticproperties, failureproperties)
V3 = Lamina(0.25, 0, elasticproperties, failureproperties)

H0 = Lamina(0.25, 0, elasticproperties, failureproperties)
H1 = Lamina(0.25, 0, elasticproperties, failureproperties)
H2 = Lamina(0.25, 0, elasticproperties, failureproperties)
H3 = Lamina(0.25, 0, elasticproperties, failureproperties)

# create the laminas list, for the laminate function:
LaminasV = [V0, V1, V2, V3]
LaminasH = [H0, H1, H2, H3]

# creating the laminate object:
LaminateH = Laminate(LaminasH)
LaminateV = Laminate(LaminasV)

TStringer = TStringer(LaminateH, LaminateV, 30, 30)
print('EI equivalent:', np.round(TStringer.FindEIEquivalent()*1e-6, 1), 'Nm^2')
print('Ybar:', TStringer.ybar)