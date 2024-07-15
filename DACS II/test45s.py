from Toolbox.Member import Member
from Toolbox.Laminate import LaminateBuilder
import numpy as np
np.set_printoptions(precision=2, linewidth=1000)
# Build a laminate
#1st: symmetrc, 2nd: copy center, number muciplicity
Laminate = LaminateBuilder([45, -45, 45, -45, 45, -45, 45, -45, 45, -45, 45, -45], True, True, 1)
print(Laminate.Ex)
# >>>>>>> 82b0b2dd12e2ccb3e857e22a392d0e4efdf24b56
# Instantiate the Member object
Member = Member(Laminate)

# Assign loading ratio and find critical load:
Member.Loads = [-1, 0, 0, 0, 0, 0]
print(np.round(Member.PanelFI(), 0))
print(Member.NormalBucklingFI())

# define impact force:
impactforce = Member.impactforce(150, 100)
deflection = Member.compute_deflection(impactforce, 150, 100, 150, 100)
print('Impact force: ', np.round(impactforce, 1), 'N')
print('Deflection at this load: ', np.round(deflection, 1), 'mm')

Member.plot_Taurz(0, 12)
Member.plot_ForceEquilibrium(0, 8)
azimuth = 0
delamination_lengths = Member.DelaminationAnalysis(azimuth, 20)
Member.plot_delamination(delamination_lengths, azimuth)

# Generate the damaged region:
damagedregion = Member.GenerateDamagedRegion(50)
print(Member.damagedregion.E_reduced)
print('Major, minor axes:', Member.Major_Minor_axes())
print(Member.CalculateCAI(0.7, 50))