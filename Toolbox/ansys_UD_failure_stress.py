import Data.MP

from Toolbox import structural_entity
from Toolbox.laminate import laminate_builder

laminate2 = laminate_builder([0], True, True, 1, 'T700')

load = 1 #N/mm
height = laminate2.h
stress = load/height
print(stress)
laminate2.Loads = [1,0,0,0,0,0]
print(laminate2.failure_analysis())

failure_stress = stress/laminate2.failure_analysis()
print('failure stress for UD T700 fiber:', failure_stress)