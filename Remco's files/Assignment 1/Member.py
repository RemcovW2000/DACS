import numpy as np
from Laminate import Laminate
from Lamina import Lamina

# The member class must do the following:

# Take the width, orientation and laminate of a member
# Calculate the:
# - Area
# - EA equivalent
# - EI equivalent around any axis? -> at least x and y
# - X and Y position

class Member(object):
    def __init__(self, Laminate, Width, Orientation):
        self.Laminate = Laminate
        self.Width = Width
        self.Orientation = Orientation


    def CalculateEI(self):
        self.Laminate.ABD_matrix_inverse
