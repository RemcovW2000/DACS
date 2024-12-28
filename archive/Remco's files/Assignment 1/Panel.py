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

class Panel(object):
    def __init__(self, Length = None, Positions = np.array([[0, 0], [0, 0]]), Laminate = None):
        self.Laminate = Laminate
        self.Position = Positions # A 2x2 numpy array, [[x y],[x, y]]
        self.Stringers = [0, 0]
        self.Length = Length
        self.EA = 0

        self.CalculateEA()

    def CalculateEA(self):          # Independent of position
        self.EA = 0  # Placeholder
        return self.EA

    def calculate_EI(self, ybar):    # give ybar as an input
        self.EI = 0  # Placeholder
        return self.EI