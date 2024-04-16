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
        # TODO: think about whether stress needs to be stored here?


# --------------------------------------------------
# Skin 1:
# --------------------------------------------------
# create laminate object
Laminate_a = LaminateBuilder([45, -45, 90, 0 , 0], True, True, 1)
# Skin 1
Skin_1 = Skin(Laminate_a)


# --------------------------------------------------
# Skin 2:
# --------------------------------------------------
# create laminate object
Laminate_b = LaminateBuilder([45, -45, 90, 0 , 0], True, True, 1)
# Skin 1
Skin_2 = Skin(Laminate_b)


# -------------------------------------------------
# Skin 3:
# -------------------------------------------------
# create laminate object
Laminate_c = LaminateBuilder([45, -45, 90, 0 , 0], True, True, 1)
# Skin 1
Skin_3 = Skin(Laminate_c)



