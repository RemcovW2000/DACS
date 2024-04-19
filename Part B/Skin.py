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
        self.A           = None
        self.Diameter    = 6000 #mm
        # TODO: think about whether stress needs to be stored here?

    def CalculateA(self):
        theta = self.stop - self.start
        if theta < 0:
            theta = self.stop + 360 - self.start
        arcl = self.Diameter * np.pi * theta/360
        self.A = self.Laminate.h * arcl
        return self.A
    
    def Calculatempl(self):
        self.CalculateA()
        self.mpl = self.rho * self.A
        return self.mpl
    
    def StaticMomentSkin(self):
        start = np.radians(self.start)
        stop  = np.radians(self.stop)
        self.static_moment = (np.cos(start) - np.cos(stop)) * self.Laminate.h * (self.Diameter/2)*2 
        return self.static_moment


# --------------------------------------------------
# Skin 1:
# --------------------------------------------------
# create laminate object

Laminate_a = LaminateBuilder([45, -45, 45, -45, 0, 45,0, -45, 45, -45, 0, 90, 0], True, True, 3)
# Skin 1
Skin_compression = Skin(Laminate_a)

# --------------------------------------------------
# Skin 2:
# --------------------------------------------------
# create laminate object
Laminate_b = LaminateBuilder([-45, 45, -45, 45, -45, 45, -45, 45 ,-45, 0, 45, -45, 0, 45, -45, 45, -45, 45, 0, 90, 0 , 90], True, True, 1)
# Skin 1
Skin_shear = Skin(Laminate_b)

# -------------------------------------------------
# Skin 3:
# -------------------------------------------------
# create laminate object
Laminate_c = LaminateBuilder([ 0, 45, 0, -45, 0, 45, -45, 0, 90, 0], True, True, 2)
# Skin 1
Skin_tension = Skin(Laminate_c)
