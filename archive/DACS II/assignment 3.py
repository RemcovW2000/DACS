import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from Toolbox import lamina
from Toolbox import laminate
from Toolbox.laminate import laminate_builder
laminate = laminate_builder([45, -45, 0, 0, 90, -45, 45], True, True, 1)
E_theta = laminate.calculate_equivalent_properties()[0][0]
print('E_theta = ', E_theta)
w = 1 # width
q1 = 1
# r = variable
t_total = 3.64
a = 7
b = a+t_total
E_r = 11.2*1000
nu_theta_r = 0.3
G_theta_r = 5*1000
d = 4.36 # mm
R = 7 # mm
S = 100 # Mpa
Yt = 70 # Mpa

gamma = np.sqrt(1 + (E_theta / E_r) * (1 - 2 * nu_theta_r) + (E_theta / G_theta_r))
q1_expr = (2 / gamma) * (1 - (a / b) ** gamma) + (1 + (a / b) ** gamma) * np.log(a / b)
k = np.sqrt(E_theta / E_r)
q2 = (1 - (a / b)**2) / 2 - (k / (k + 1)) * ((1 - (a / b)**(k + 1))**2 / (1 - (a / b)**(2 * k))) + \
          (k * (a / b)**2 / (k - 1)) * ((1 - (a / b)**(k - 1))**2 / (1 - (a / b)**(2 * k)))


def Tau(r, theta, F):
    Tau = (F / (r * w * q1)) * (((r / b) ** gamma) + ((a / b) ** gamma) * ((b / r) ** gamma) - 1 - ((a / b) ** gamma)) * np.sin(theta)
    return Tau

def sigma_rs(r, theta, F):
    sigma_r = (F * np.cos(theta) / (r * w * q1)) * (((r ** gamma) / (b ** gamma)) + ((a ** gamma) / (r ** gamma)) - 1 - (a ** gamma / b ** gamma))
    return sigma_r

def sigma_rm(r, theta, F):
    M = (d + np.cos(theta) * (R + t_total/2))*F
    sigma_rm = -(M / (b ** 2 * w * q2)) * (
            1 -
            (1 - (a / b) ** (k + 1)) / (1 - (a / b) ** (2 * k)) * (r / b) ** (k - 1) -
            (1 - (a / b) ** (k - 1)) / (1 - (a / b) ** (2 * k)) * (a / b) ** (k + 1) * (b / r) ** (k + 1)
    )
    return sigma_rm
def Failure_criterion(r, theta, F):
    sigma_r = sigma_rm(r, theta, F) + sigma_rs(r, theta, F)
    Tau1 = Tau(r, theta, F)

    FI = (Tau1 / S)**2 + (sigma_r/Yt)**2
    return FI


def find_maximum_FI(r_values, theta_values ,F):
    # Create meshgrid for r and theta
    R, Theta = np.meshgrid(r_values, theta_values)

    # Calculate Failure Criterion values for the meshgrid
    FI_values = Failure_criterion(R, Theta, F)

    # Find the index of the maximum FI value
    max_index = np.unravel_index(np.argmax(FI_values, axis=None), FI_values.shape)

    # Get the r and theta values corresponding to the maximum FI value
    r_max = R[max_index]
    theta_max = Theta[max_index]
    max_FI = FI_values[max_index]

    return r_max, theta_max, max_FI


def find_F_for_max_FI_equal_one(r_values, theta_values):
    def objective(F):
        _, _, max_FI = find_maximum_FI(r_values, theta_values, F)
        return np.abs(max_FI - 1)  # We want max_FI to be 1

    result = minimize_scalar(objective, bounds=(50, 200), method='bounded')
    return result.x


# Range of r values
r_values = np.arange(7.26, 10.38 + 0.26, 0.26)  # From 7.26 to 10.38 in intervals of 0.26

# Generate theta values
theta_values = np.linspace(0, 0.5 * np.pi, 100)


# Find the value of F for which the maximum FI is 1
optimal_F = find_F_for_max_FI_equal_one(r_values, theta_values)
print(f"The value of F for which the maximum FI is 1 is: {optimal_F:.4f}")
F = optimal_F*1.3
maxdata = find_maximum_FI(r_values, theta_values, F)
print('max FI = ',maxdata[2],' is found at r =', maxdata[0], ' theta = ', maxdata[1])

# Create meshgrid for r and theta
R1, Theta = np.meshgrid(r_values, theta_values)

# Calculate Failure Criterion values for the meshgrid
FI_values = Failure_criterion(R1, Theta, F)


# Plotting
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Create the surface plot
surf = ax.plot_surface(R1, Theta, FI_values, cmap='viridis')

# Add color bar which maps values to colors
fig.colorbar(surf, shrink=0.5, aspect=5)

# Add a translucent horizontal plane at z = 1
x_plane = np.linspace(r_values.min(), r_values.max(), 100)
y_plane = np.linspace(theta_values.min(), theta_values.max(), 100)
X_plane, Y_plane = np.meshgrid(x_plane, y_plane)
Z_plane = np.ones_like(X_plane)
ax.plot_surface(X_plane, Y_plane, Z_plane, color='blue', alpha=0.3)


# Change the view angle
ax.view_init(elev=0, azim=90)  # Adjust the elevation and azimuthal angle as needed


# Labels and title
ax.set_xlabel('r')
ax.set_ylabel('Theta (radians)')
ax.set_zlabel('Failure Criterion')
ax.set_title('Failure Criterion as a function of r and Theta, F={} N/mm'.format(np.round(F, 2)))
newlist = []
list = [0.4, 0.75, 0.82, 0.85, 0.8, 0.7, 0.42]
r_values = np.arange(7.78, 19.34 + 0.26, 0.26)
for index, angle in enumerate(list):
    distance = angle*r_values[index]
    newlist.append(distance)

print(newlist)
plt.show()