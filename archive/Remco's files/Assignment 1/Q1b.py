import numpy as np
from Laminate import Laminate
from Lamina import Lamina
import matplotlib.pyplot as plt

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

# Defining lamina:
s0 = Lamina(0.125, 0, elasticproperties, failureproperties)
s1 = Lamina(0.125, 0, elasticproperties, failureproperties)
s2 = Lamina(0.125, 90, elasticproperties, failureproperties)
s3 = Lamina(0.125, 30, elasticproperties, failureproperties)
s4 = Lamina(0.125, 90, elasticproperties, failureproperties)

laminas = [s0, s1, s2, s3, s4]

# Creating laminte object
laminate = Laminate(laminas)

# Defining loads
laminate.Loads = np.array([[0.02], [18], [0], [18000], [0], [0]])

# First we have to calculate strains:
laminate.calculate_lamina_strains()

# Now we print the stresses and strains
print('S0: ',s0.calculate_principal_sigma_epsilon())
print('S1: ',s1.calculate_principal_sigma_epsilon())
print('S2: ', s2.calculate_principal_sigma_epsilon())
print('S3: ', s3.calculate_principal_sigma_epsilon())
print('S4: ', s4.calculate_principal_sigma_epsilon())