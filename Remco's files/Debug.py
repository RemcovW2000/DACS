import numpy as np
from Laminate import Laminate
from Lamina import Lamina
import matplotlib.pyplot as plt
from tqdm import tqdm

np.set_printoptions(linewidth=300, precision=2)

# now we test the code:
E1 = 145300     # From assignment
E2 = 8500       # From assignment
G12 = 4580      # From assignment
v12 = 0.31      # From assignment

elasticproperties = [E1, E2, G12, v12]

E11f = 500000   # TBD
v21f = 0.1      # TBD
msf = 1.1       # TBD
R11t = 1932     # From assignment
R11c = 1480     # From assignment
yt = 108        # From assignment
yc = 220        # From assignment
S = 132.8       # From assignment

failureproperties = [E11f, v21f, msf, R11t, R11c, yt, yc, S]
s0 = Lamina(0.25, 0, elasticproperties, failureproperties)
s1 = Lamina(0.25, 90, elasticproperties, failureproperties)
s2 = Lamina(0.25, 45, elasticproperties, failureproperties)
s3 = Lamina(0.25, -45, elasticproperties, failureproperties)
s4 = Lamina(0.25, -45, elasticproperties, failureproperties)
s5 = Lamina(0.25, 45, elasticproperties, failureproperties)
s6 = Lamina(0.25, 90, elasticproperties, failureproperties)
s7 = Lamina(0.25, 0, elasticproperties, failureproperties)
s8 = Lamina(0.25, 0, elasticproperties, failureproperties)
s9 = Lamina(0.25, 90, elasticproperties, failureproperties)
s10 = Lamina(0.25, 45, elasticproperties, failureproperties)
s11 = Lamina(0.25, -45, elasticproperties, failureproperties)
s12 = Lamina(0.25, -45, elasticproperties, failureproperties)
s13 = Lamina(0.25, 45, elasticproperties, failureproperties)
s14 = Lamina(0.25, 90, elasticproperties, failureproperties)
s15 = Lamina(0.25, 0, elasticproperties, failureproperties)

# create the laminas list, for the laminate function:
laminas = [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15]

# creating the laminate object:
laminate = Laminate(laminas)


E22vsE12FPF, E22vsE12LPF, S22vsS12FPF, S22vsS12LPF = laminate.ProduceFailureEnvelope(0.1)

# Unzip the list of tuples into two lists, x and y
Fx, Fy = zip(*S22vsS12FPF)
Lx, Ly = zip(*S22vsS12LPF)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(Lx, Ly, color='red')
plt.plot(Fx, Fy, color='blue')

plt.title('FPF and LPF stresses')
plt.xlabel('X values')
plt.ylabel('Y values')
plt.grid(True)
plt.axis('equal')
plt.show()