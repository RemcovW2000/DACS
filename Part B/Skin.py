import numpy as np
from Laminate import Laminate
from Lamina import Lamina
import MP

class Skin:
    def __init__(self, Laminate=None):
        self.Laminate    = Laminate
        self.start       = 0
        self.stop        = 0
        #TODO: think about whether stress needs to be stored here?

# --------------------------------------------------
# Skin 1:
# --------------------------------------------------
a0 = Lamina(MP.t, 45, MP.elasticproperties, MP.failureproperties)
a1 = Lamina(MP.t, -45, MP.elasticproperties, MP.failureproperties)
a2 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
a3 = Lamina(MP.t, 90, MP.elasticproperties, MP.failureproperties)
a4 = Lamina(MP.t, 90, MP.elasticproperties, MP.failureproperties)
a5 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
a6 = Lamina(MP.t, -45, MP.elasticproperties, MP.failureproperties)
a7 = Lamina(MP.t, 45, MP.elasticproperties, MP.failureproperties)

# create lamina lists
Laminas_a = [a0]
# create laminate object
Laminate_a = Laminate(Laminas_a)
# Skin 1
Skin_a = Skin(Laminate_a)

# --------------------------------------------------
# Skin 2:
# --------------------------------------------------
b0 = Lamina(MP.t, 45, MP.elasticproperties, MP.failureproperties)
b1 = Lamina(MP.t, -45, MP.elasticproperties, MP.failureproperties)
b2 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
b3 = Lamina(MP.t, 90, MP.elasticproperties, MP.failureproperties)
b4 = Lamina(MP.t, 90, MP.elasticproperties, MP.failureproperties)
b5 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
b6 = Lamina(MP.t, -45, MP.elasticproperties, MP.failureproperties)
b7 = Lamina(MP.t, 45, MP.elasticproperties, MP.failureproperties)


# create lamina lists
Laminas_b = [b0]
# create laminate object
Laminate_b = Laminate(Laminas_b)
# Skin 1
Skin_b = Skin(Laminate_b)

# -------------------------------------------------
# Skin 3:
# -------------------------------------------------
c0 = Lamina(MP.t, 45, MP.elasticproperties, MP.failureproperties)
c1 = Lamina(MP.t, -45, MP.elasticproperties, MP.failureproperties)
c2 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
c3 = Lamina(MP.t, 90, MP.elasticproperties, MP.failureproperties)
c4 = Lamina(MP.t, 90, MP.elasticproperties, MP.failureproperties)
c5 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
c6 = Lamina(MP.t, -45, MP.elasticproperties, MP.failureproperties)
c7 = Lamina(MP.t, 45, MP.elasticproperties, MP.failureproperties)

# create lamina lists
Laminas_c = [c0]
# create laminate object
Laminate_c = Laminate(Laminas_c)
# Skin 1
Skin_c = Skin(Laminate_c)



