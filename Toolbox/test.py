import numpy as np

v1 = np.array([2, 2])
v2 = np.array([3, 3])

v3 = v2 + v1
norm = v3 / np.linalg.norm(v2+v1)
print(np.sqrt(norm[0]**2 + norm[1]**2))
print(v3)
print(norm)