import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
Rt = 1000e6
Rc = 800e6
v21 = 0
v21f = 0
msigmaf = 1.2
E11 = 140e9
E11f = 150e9
factor = v21 - v21f * msigmaf * (E11 / E11f)

# Range of S2 and S3 values
S2 = np.linspace(-100e9, 10, 400)
S3 = np.linspace(-10, 10, 400)

# Meshgrid for S2 and S3
S2, S3 = np.meshgrid(S2, S3)

# Solving for S1 when f_EFF = 1 for both conditions
# For S1 > 0 (using Rt)
S1_positive = Rt * (1 + factor * (S2 + S3))

# For S1 < 0 (using Rc)
S1_negative = Rc * (1 + factor * (S2 + S3))

# Plotting
fig = plt.figure(figsize=(14, 7))

# S1 > 0
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(S2, S3, S1_positive, color='blue', alpha=0.5, rstride=100, cstride=100)
ax1.set_title('S1 > 0')
ax1.set_xlabel('S2')
ax1.set_ylabel('S3')
ax1.set_zlabel('S1')
ax1.view_init(30, 120)

# S1 < 0
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(S2, S3, S1_negative, color='red', alpha=0.5, rstride=100, cstride=100)
ax2.set_title('S1 < 0')
ax2.set_xlabel('S2')
ax2.set_ylabel('S3')
ax2.set_zlabel('S1')
ax2.view_init(30, 120)

plt.tight_layout()
plt.show()