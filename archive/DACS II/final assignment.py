from Toolbox.member import Member
from Toolbox.laminate import laminate_builder
import numpy as np
np.set_printoptions(precision=2, linewidth=1000)
# Build a laminate
#1st: symmetrc, 2nd: copy center, number muciplicity
Laminate = laminate_builder([45, -45, 45, -45, 45, -45, 45, -45, 0, 90, 0], True, True, 1)
print(Laminate.Ex)
# >>>>>>> 82b0b2dd12e2ccb3e857e22a392d0e4efdf24b56
# Instantiate the Member object
Member = Member(Laminate)

# Assign loading ratio and find critical load:
Member.Loads = [-1, 0, 0, 0, 0, 0]
print(np.round(Member.panel_FI(), 0))
print(Member.normal_load_buckling_FI())

# define impact force:
impactforce = Member.calculate_impactForce(50, 50)
deflection = Member.compute_deflection(impactforce, 50, 50, 50, 50)
print('Impact force: ', np.round(impactforce, 1), 'N')
print('Deflection at this load: ', np.round(deflection, 1), 'mm')

Member.plot_Taurz(0, 12)
Member.plot_force_equilibrium(0, 8)
azimuth = 0
delamination_lengths = Member.delamination_analysis(azimuth, 20)
Member.plot_delamination(delamination_lengths, azimuth)

# Generate the damaged region:
damagedregion = Member.generate_damaged_region(50)
print(Member.damagedregion.E_reduced)
print('Major, minor axes:', Member.major_minor_axes())
print(Member.calculate_CAI(0.7, 50))
