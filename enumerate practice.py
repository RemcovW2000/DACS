import numpy as np
S1 = 100e6
S2 = 0
S3 = 0

E11 = 140e9     # e modulus of lamina in 11 direction
E11f = 100      # Fiber e modulus
msigmaf = 1.1   # Intended to capture the differences in the transverse stresses in the fiber and matrx.
v21 = 0         # possion ratio lamina
v21f = 100      # possion ratio fiber

Rt = 1000e9
Rc = 800e9

if S1 > 0:
    f_EFF = (1 / Rt) * (S1 - (v21 - v21f * msigmaf * (E11 / E11f)) * (S2 + S3))
elif S1 < 0:
    f_EFF = (1 / Rc) * (S1 - (v21 - v21f * msigmaf * (E11 / E11f)) * (S2 + S3))

array = np.array([[1], [2], [3]])

print(array[1])