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
R11c = 1932     # TBD
yt = 108        # From assignment
yc = 108        # TBD
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

# Making the list of angles:
angles = np.linspace(1, 360, 5)

Failure_shear_strains = np.array([])
Failure_transverse_strains = np.array([])


for angle in tqdm(angles):
    # For the angle, we create the loading ratio in which direction we will increment the load:
    loadingratio = np.array([[0],
                             [np.sin(np.deg2rad(angle))],
                             [np.cos(np.deg2rad(angle))],
                             [0],
                             [0],
                             [0]])

    # Then we calulate the FPF and LPF load:
    FPF_load, LPF_load, FailureLoads, FailureStrains = laminate.ProgressiveDamageAnalysis(loadingratio, 1)

    # Finally we have to reset the failurestate of the laminate:
    laminate.ResetFailureState()

    # Now we need to plot the global failure strains, we'll plot the shear and transverse strain (row 1 and 2)
    Failure_shear_strains = np.append(Failure_shear_strains, FailureStrains[1])
    Failure_transverse_strains = np.append(Failure_transverse_strains, FailureStrains[2])

laminate.PlotFailureEnvelope()
# Create a scatter plot of the failure strains:
plt.scatter(Failure_transverse_strains, Failure_shear_strains, color='red', label='strains')

# Optionally, you can add titles and labels
plt.title('Failure strains')
plt.xlabel('shears')
plt.ylabel('transverses')

# Show the plot
plt.show()
