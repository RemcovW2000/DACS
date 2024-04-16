import numpy as np
from Laminate import Laminate
from Laminate import LaminateBuilder
from Lamina import Lamina
import MP

class Skin:
    def __init__(self, Laminate=None):
        self.Laminate    = Laminate
        self.start       = 0
        self.stop        = 0
        self.mpl         = 0
        self.rho         = 0.00161
        self.A           = 0
        self.Diameter    = 0
        # TODO: think about whether stress needs to be stored here?
        self.CalculateA()
    def CalculateA(self):
        arcl = self.Diameter * np.pi * (self.stop-self.start)/360
        self.A = self.Laminate.h * arcl
        return self.A
    def Calculatempl(self):
        self.CalculateA()
        self.mpl = self.rho * self.A
        return self.mpl


# --------------------------------------------------
# Skin 1:
# --------------------------------------------------
# create laminate object
Laminate_a = LaminateBuilder([45, -45, 90, 0 , 0], True, True, 1)
# Skin 1
Skin_compression = Skin(Laminate_a)


# --------------------------------------------------
# Skin 2:
# --------------------------------------------------
# create laminate object
Laminate_b = LaminateBuilder([45, -45, 45, -45, 45, -45, 45, -45, 0], True, True, 5)
# Skin 1
Skin_shear = Skin(Laminate_b)


# -------------------------------------------------
# Skin 3:
# -------------------------------------------------
# create laminate object
Laminate_c = LaminateBuilder([45, -45, 90, 0 , 0], True, True, 1)
# Skin 1
Skin_tension = Skin(Laminate_c)



