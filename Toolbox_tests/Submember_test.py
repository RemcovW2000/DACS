from Data.Airfoils import airfoilcoords
from scipy.interpolate import interp1d
from Toolbox.laminate import laminate_builder
import matplotlib.pyplot as plt
from Toolbox.damaged_region import *
from Toolbox.member import Member
from scipy.integrate import simps
from matplotlib.colors import Normalize, TwoSlopeNorm
from scipy.linalg import lstsq
import math
from Toolbox.section import section
from Toolbox.helperclasses import boom, segment

mainlaminate = laminate_builder([0, 90, 45, -45], True, True, 1)

reinforcedlaminate = laminate_builder([0, 90, -45, 45, 0, 0, 0], True, True, 1)

member = Member(mainlaminate, reinforcedlaminate)
member.subpanel = reinforcedlaminate
member.subpanel_start = 0
member.subpanel_end = 20

print(member.get_Ex(30))

import numpy as np

# Example linspace and extra value
x_values = np.linspace(0, 10, 11)  # [0, 1, 2, ..., 10]
x = 4.5  # Value to insert

# Find the index where x should be inserted
index = np.searchsorted(x_values, x)

# Insert x into the array
new_x_values = np.insert(x_values, index, x)

print("Original array:", x_values)
print("New array:", new_x_values)