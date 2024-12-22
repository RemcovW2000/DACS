import numpy as np
from Laminate import Laminate
from Lamina import Lamina
import matplotlib.pyplot as plt


# Assigning elastic properties:
E1 = 145.3e3     # From assignment
E2 = 8.5e3       # From assignment
G12 = 4.58e3      # From assignment
v12 = 0.31      # From assignment

elasticproperties = [E1, E2, G12, v12]

# properties needed for failure analysis
E11f = 230e3
v21f = 0.2
msf = 1.1
R11t = 1932     # From assignment
R11c = 1480     # From assignment
yt = 108        # From assignment
yc = 220        # From assignment
S = 132.8       # From assignment

failureproperties = [E11f, v21f, msf, R11t, R11c, yt, yc, S]

# Defining the lamina:
s0 = Lamina(0.125, 15, elasticproperties, failureproperties)
s1 = Lamina(0.125, 30, elasticproperties, failureproperties)
s2 = Lamina(0.125, -30, elasticproperties, failureproperties)
s3 = Lamina(0.125, 75, elasticproperties, failureproperties)
s4 = Lamina(0.125, 75, elasticproperties, failureproperties)
s5 = Lamina(0.125, 15, elasticproperties, failureproperties)
s6 = Lamina(0.125, 30, elasticproperties, failureproperties)
s7 = Lamina(0.125, -30, elasticproperties, failureproperties)
s8 = Lamina(0.125, 75, elasticproperties, failureproperties)
s9 = Lamina(0.125, 75, elasticproperties, failureproperties)
s10 = Lamina(0.125, 75, elasticproperties, failureproperties)
s11 = Lamina(0.125, 75, elasticproperties, failureproperties)
s12 = Lamina(0.125, -30, elasticproperties, failureproperties)
s13 = Lamina(0.125, 30, elasticproperties, failureproperties)
s14 = Lamina(0.125, 15, elasticproperties, failureproperties)
s15 = Lamina(0.125, 75, elasticproperties, failureproperties)
s16 = Lamina(0.125, 75, elasticproperties, failureproperties)
s17 = Lamina(0.125, -30, elasticproperties, failureproperties)
s18 = Lamina(0.125, 30, elasticproperties, failureproperties)
s19 = Lamina(0.125, 15, elasticproperties, failureproperties)

# Producing the lamina list
laminas = [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18, s19]

# creating the object
laminate = Laminate(laminas)

# Finding the engineering constants:
print(laminate.CalculateEquivalentProperties()[1])