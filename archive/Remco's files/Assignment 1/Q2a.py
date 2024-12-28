import numpy as np
from Laminate import Laminate
from Lamina import Lamina
import matplotlib.pyplot as plt
from tqdm import tqdm

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
s0 = Lamina(0.125, 0, elasticproperties, failureproperties)
s1 = Lamina(0.125, 90, elasticproperties, failureproperties)
s2 = Lamina(0.125, 45, elasticproperties, failureproperties)
s3 = Lamina(0.125, -45, elasticproperties, failureproperties)
s4 = Lamina(0.125, -45, elasticproperties, failureproperties)
s5 = Lamina(0.125, 45, elasticproperties, failureproperties)
s6 = Lamina(0.125, 90, elasticproperties, failureproperties)
s7 = Lamina(0.125, 0, elasticproperties, failureproperties)
s8 = Lamina(0.125, 0, elasticproperties, failureproperties)
s9 = Lamina(0.125, 90, elasticproperties, failureproperties)
s10 = Lamina(0.125, 45, elasticproperties, failureproperties)
s11 = Lamina(0.125, -45, elasticproperties, failureproperties)
s12 = Lamina(0.125, -45, elasticproperties, failureproperties)
s13 = Lamina(0.125, 45, elasticproperties, failureproperties)
s14 = Lamina(0.125, 90, elasticproperties, failureproperties)
s15 = Lamina(0.125, 0, elasticproperties, failureproperties)

# create the laminas list, for the laminate function:
laminas = [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15]
# laminas = [s0, s1, s2,s3, s4, s5, s6, s7]
# laminas = [s0, s1, s6, s7]
# creating the laminate object:
laminate = Laminate(laminas)


E22vsE12FPF, E22vsE12LPF, S22vsS12FPF, S22vsS12LPF, Strains = laminate.produce_failure_envelope(0.1)
# Unzip the list of tuples into two lists, x and y
Fx, Fy = zip(*S22vsS12FPF)
Lx, Ly = zip(*S22vsS12LPF)

# Flatten the list of tuples
flat_data = [point for sublist in Strains for point in sublist]
E22 = []
E12 = []
for transverse, shear in flat_data:
    E22.append(transverse)
    E12.append(shear)



# Plotting
# fig, axs = plt.subplots(2, figsize=(6, 8), gridspec_kw={'hspace': 0.3})
#
# axs[0].plot(Lx, Ly, color='red', label = 'Last ply failure')
# axs[0].plot(Fx, Fy, color='blue', label = 'First ply damage')
# axs[0].set_title('First ply failure, last ply failure, stresses')
# axs[0].set_xlabel('Ny (N/mm)')
# axs[0].set_ylabel('Ns (N/mm)')
# axs[0].legend()
# axs[0].grid(True)
# axs[0].axis('equal')

# axs[1].scatter(E22, E12, color='red', label = 'Failure point', s=0.1, marker='.')
# axs[1].set_title('All failures, strains')
# axs[1].set_xlabel('Ey')
# axs[1].set_ylabel('Ex')
# axs[1].legend()
# axs[1].axis('equal')
# axs[1].grid(True)

plt.plot(Lx, Ly, color='red', label = 'Last ply failure')
plt.plot(Fx, Fy, color='blue', label = 'First ply damage')
plt.xlabel('Ny (N/mm')
plt.ylabel('Ns (N/mm)')
plt.axis('equal')
plt.title('Stresses at which failure occurs, max stress criteria')
plt.legend()
plt.grid(True)
plt.savefig('Q2a stresses, FPF, LPF, max_stress', dpi=1200)
plt.show()