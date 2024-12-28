import sys
import os
print("Current working directory:", os.getcwd())
print("Script directory:", os.path.dirname(os.path.abspath(__file__)))

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
print("Updated sys.path:", sys.path)




from Toolbox.member import Member
from Toolbox.laminate import laminate_builder
import numpy as np
np.set_printoptions(precision=2, linewidth=1000)
# Build a laminate
#1st: symmetrc, 2nd: copy center, number muciplicity
Laminate = laminate_builder([45, -45, 45, -45, 45, -45, 45, -45, 0, 90, 0], True, True, 1)
# >>>>>>> 82b0b2dd12e2ccb3e857e22a392d0e4efdf24b56
# Instantiate the Member object
Member = Member(Laminate)

# Assign loading ratio and find critical load:
Member.Loads = [-1, 0, 0, 0, 0, 0]
print(np.round(Member.panel_FI(), 0))
print(Member.normal_load_buckling_FI())

# define impact force:
impactforce = Member.calculate_impactForce(150, 100)
deflection = Member.compute_deflection(impactforce, 150, 100, 150, 100)
print('Impact force: ', np.round(impactforce, 1), 'N')
print('Deflection at this load: ', np.round(deflection, 1), 'mm')

Member.plot_Taurz(0, 12)
Member.plot_force_equilibrium(0, 8)
azimuth = 0
delamination_lengths = Member.delamination_analysis(azimuth, 20)
print(max(delamination_lengths))
Member.plot_delamination(delamination_lengths, azimuth)

# Generate the damaged region:
damagedregion = Member.generate_damaged_region(50)
print('Major, minor axes:', Member.major_minor_axes())
print(Member.calculate_CAI(0.7, 50))
