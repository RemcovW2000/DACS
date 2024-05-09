import numpy as np
from tqdm import tqdm
import copy
from Laminate import Laminate
from Lamina import Lamina
import MP
from Laminate import LaminateBuilder



laminateshear = LaminateBuilder([0, 0, 0, 0], False, False, 1)
print(laminateshear.D_matrix)

s0 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
s1 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
s2 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
s3 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)

laminas = [s0, s1, s2, s3]
laminate = Laminate(laminas)
