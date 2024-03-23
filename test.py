import numpy as np
import matplotlib.pyplot as plt

# Define the range for theta
theta_deg = np.linspace(0, 180, 360)
theta_rad = np.radians(theta_deg)

# Material properties
E1 = 140E9
G12 = 12E9
v12 = 0.3
E2 = 8E9

# Preallocate arrays for storing computed values
Ex = np.zeros_like(theta_rad)
Ey = np.zeros_like(theta_rad)
vxy = np.zeros_like(theta_rad)
Gxy = np.zeros_like(theta_rad)
nxs = np.zeros_like(theta_rad)
nys = np.zeros_like(theta_rad)

# Compute the values for each theta
for i, theta in enumerate(theta_rad):
    m = np.cos(theta)
    n = np.sin(theta)

    Ex[i] = 1 / ((m ** 4 / E1) + ((1 / G12) - 2 * (v12 / E1)) * m ** 2 * n ** 2 + ((n ** 4) / E2))
    vxy[i] = Ex[i] * (v12 / E1 * (m ** 4 + n ** 4) - (1 / E1 + 1 / E2 - 1 / G12) * m ** 2 * n ** 2)
    Ey[i] = 1 / (n ** 4 / E1 + (1 / G12 - 2 * v12 / E1) * m ** 2 * n ** 2 + m ** 4 / E2)
    Gxy[i] = 10*(1 / (2 * (2 / E1 + 2 / E2 + 4 * v12 / E1 - 1 / G12) * m ** 2 * n ** 2 + 1 / G12 * (m ** 4 + n ** 4)))
    nxs[i] = Ex[i] * ((2 / E1 + 2 * v12 / E1 - 1 / G12) * m ** 3 * n - (2 / E2 + 2 * v12 / E1 - 1 / G12) * m * n ** 3)
    nys[i] = Ey[i] * ((2 / E1 + 2 * v12 / E1 - 1 / G12) * m * n ** 3 - (2 / E2 + 2 * v12 / E1 - 1 / G12) * m ** 3 * n)

# Plotting
plt.figure(figsize=(20, 10))

plt.subplot(2, 2, 1)
plt.plot(theta_deg, Ex, label='Ex')
plt.plot(theta_deg, Ey, label='Ey')
plt.plot(theta_deg, Gxy, label = 'Gxy', color='green')
plt.title('Ex and Ey vs Theta')
plt.xlabel('Theta (degrees)')
plt.ylabel('Elastic Modulus (Pa)')
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(theta_deg, vxy, color='orange')
plt.title('vxy vs Theta')
plt.xlabel('Theta (degrees)')
plt.ylabel('Poisson\'s Ratio')

plt.subplot(2, 2, 3)
plt.plot(theta_deg, nxs, label='nxs')
plt.plot(theta_deg, nys, label='nys')
plt.title('nxs and nys vs Theta')
plt.xlabel('Theta (degrees)')
plt.ylabel('Shear Coupling Coefficient')

plt.tight_layout()
plt.show()

