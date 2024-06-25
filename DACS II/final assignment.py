from Toolbox.Member import Member
from Toolbox.Laminate import LaminateBuilder
import numpy as np
np.set_printoptions(precision=2, linewidth=1000)
# Build a laminate
#1st: symmetrc, 2nd: copy center, number muciplicity 
Laminate = LaminateBuilder([45, -45, 0, 90], True, True, 2)
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

# delamination_lengths = Member.DelaminationAnalysis(0, 20)
# Member.plot_delamination(delamination_lengths)

# Generate the damaged region:
damagedregion = Member.GenerateDamagedRegion()
print(np.round(damagedregion.E_reduced, 2))
print(np.round(Member.panel.Ex, 2))
Member.Major_Minor_axes()
print(Member.CalculateCAI(0.7))