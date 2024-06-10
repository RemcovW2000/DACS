import numpy as np
import copy
from Toolbox import MP
from Toolbox.Laminate import LaminateBuilder
from Toolbox.Sandwich import Sandwich, Core, Member

laminate1 = LaminateBuilder([45, -45, 45, -45], True, False, 1)
laminate2 = copy.deepcopy(laminate1)

thickness = 10 #mm
core = Core(thickness, MP.HRH128props)

shearforce = 55000 #N
a = 500 #mm
b = 900 #mm
Loads = [0, 0, 0, shearforce/b, 0, 0]
# Load is around 61 N/mm
sandwich = Sandwich(laminate1, laminate2, core, Loads)

Member = Member(sandwich, Loads, a, b)
print(np.round(Member.ShearBucklingFF(), 2), 'N/mm')
np.set_printoptions(precision=1, suppress=False)
print(sandwich.SandwichWrinkingShear())
print('ABD Matrix of sandwich:')
print(sandwich.ABD_matrix)