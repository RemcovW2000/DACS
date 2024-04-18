# import numpy as np
#
# E = 69000
# m = 1
# n = 2
# r = 3000
# t = 1.754
# L = 1000
# v = 0.29
#
# Z = (L**2/(r*t))*np.sqrt(1-v**2)
# beta = n*L/(m*np.pi*r)
# D = E*t**3/(12*(1-v**2))
# phi = np.sqrt(r/t)/16
# gamma = 1-0.731*(1-np.e**(-phi))
# kx = m**2*(1+beta**2)**2 + 12*(gamma*Z)**2/(np.pi**4* m**2*(1+beta**2)**2)
# Mc = kx * (np.pi ** 3 * D * r ** 2) / L ** 2
#
# print(Mc * 1e-9)
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Define parameters
E = 69000
r = 3000
t = 10
L = 1000
v = 0.3

# Define functions
def Z(L, r, t, v):
    return (L**2 / (r * t)) * np.sqrt(1 - v**2)

def beta(n, L, m, r):
    return n * L / (m * np.pi * r)

def D(E, t, v):
    return E * t**3 / (12 * (1 - v**2))

def phi(r, t):
    return np.sqrt(r / t) / 16

def gamma(phi):
    return 1 - 0.731 * (1 - np.e**(-phi))

def kx(m, n, L, r, t, v):
    Z_val = Z(L, r, t, v)
    beta_val = beta(n, L, m, r)
    phi_val = phi(r, t)
    gamma_val = gamma(phi_val)
    return m**2 * (1 + beta_val**2)**2 + 12 * (gamma_val * Z_val)**2 / (np.pi**4 * m**2 * (1 + beta_val**2)**2)


def Mc(m, n, L, r, t, v):
    D_val = D(E, t, v)
    kx_val = kx(m, n, L, r, t, v)
    return kx_val * (np.pi**3 * D_val * r**2) / L**2

def dMc_dm(m, n, L, r, t, v):
    h = 1e-5
    return (Mc(m + h, n, L, r, t, v) - Mc(m, n, L, r, t, v)) / h

def dMc_dn(m, n, L, r, t, v):
    h = 1e-5
    return (Mc(m, n + h, L, r, t, v) - Mc(m, n, L, r, t, v)) / h

# Define the objective function to minimize
# Define the combined objective function
def combined_objective(x):
    m, n = x
    return np.abs(dMc_dm(m, n, L, r, t, v)) + np.abs(dMc_dn(m, n, L, r, t, v))

initial_guess = [2, 2]
# Perform optimization
result = minimize(combined_objective, initial_guess, bounds=((1, None), (2, None)))

# Extract optimal values
optimal_m, optimal_n = result.x

print("Optimal m:", optimal_m)
print("Optimal n:", optimal_n)
print("Minimum combined derivative:", result.fun)
# Create grid of m and n values
m_values = np.linspace(1, 6, 100)
n_values = np.linspace(2, 6, 100)
m_grid, n_grid = np.meshgrid(m_values, n_values)

# Calculate Mc for each combination of m and n
Mc_values = Mc(m_grid, n_grid, L, r, t, v)

# Plot 3D surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(m_grid, n_grid, Mc_values, cmap='viridis')
ax.set_xlabel('m')
ax.set_ylabel('n')
ax.set_zlabel('Mc')
ax.set_title('Mc as a function of m and n')
plt.show()
