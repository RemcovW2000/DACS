import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from tensorflow.python.debug.cli.command_parser import parse_time_interval


# Load the CSV file
file_path = "output_01_incomplete.csv"  # Adjust path as needed
df = pd.read_csv(file_path, skiprows=3, delimiter=",", engine="python")

# Rename columns
df.columns = ["Name", "P1", "P2", "P6", "P13","P14","P15","P16","P7", "P8", "P9", "P10"]
df = df.apply(pd.to_numeric, errors='coerce').dropna()

# Define custom function to compute the angle
def compute_relative_angle(row):
    vertical_dist_inner = row["P13"] * 2
    angle_inner = np.arctan((row["P8"] - row["P7"]) / vertical_dist_inner)

    vertical_dist_outer = row["P16"] * 2
    angle_outer = np.arctan((row["P9"] - row["P10"]) / vertical_dist_outer)
    return np.degrees(angle_inner - angle_outer)  # Ensure the returned value is a scalar

def compute_overlapping_length(row):
    overlapping_length = 150 - row['P6']
    return overlapping_length

half_span = 1.5 # meters
weight = 12 # kg
g = 9.81 # m/s^2

def compute_loadcase(n):
    force = (weight * g * 0.5) * n
    arm =  4 * half_span / (3 * np.pi)
    arm_in_sim = 0.1 # m
    moment = (arm - arm_in_sim) * force
    return moment, -force

def compute_n(row):
    force = row["P1"]
    n = abs(force/(weight * g * 0.5))
    return n

# Compute and add the new column
df["Angle"] = df.apply(compute_relative_angle, axis=1)
df["overlapping_length"] = df.apply(compute_overlapping_length, axis=1)
df["n"] = df.apply(compute_n, axis=1)

# Extract input variables and the new Angle column
input_vars = df[["n", "P6", "P14"]].values.astype(np.float64)


angle_values = df["Angle"].values


# Create interpolation function for Angle
interp_Angle = LinearNDInterpolator(input_vars, angle_values)

overlapping_length_values = 150 - df["P6"].unique()
radius_values = df["P14"].unique()

# Generate data
new_data = []
n = 2
for overlapping_length in overlapping_length_values:
    for radius in radius_values:
        p6 = 150 - overlapping_length
        angle = float(interp_Angle([n, p6, radius]))  # Interpolate
        new_data.append({"n": n, "overlapping_length": overlapping_length, "Angle": angle, "radius": radius})

# Create new DataFrame
df_loadcase = pd.DataFrame(new_data)
df_n1 = df_loadcase[df_loadcase["n"] == 1]

overlapping_length_vals_n1 = df_n1["overlapping_length"].values
radius_vals_n1 = df_n1["radius"].values
angle_vals_n1 = df_n1["Angle"].values

radius_vals = df_loadcase["radius"].unique()
overlapping_length_vals = df_loadcase["overlapping_length"].unique()

# Create a 3D surface plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Create a meshgrid where radius_vals are repeated in rows, and overlapping_length_vals are repeated in columns
X, Y = np.meshgrid(radius_vals, overlapping_length_vals)
print(X, Y)
print(list(zip(np.ravel(X), np.ravel(Y))))

# Interpolate the angle values for the surface plot
Z = np.array([interp_Angle([n, 150 - overlap, radius]) for radius, overlap in zip(np.ravel(X), np.ravel(Y))])
Z = Z.reshape(X.shape)
print(interp_Angle([2, 12, 160]))
# Plot the surface
surf = ax.plot_surface(X, Y, Z, cmap='viridis')

# Labels and title
ax.set_xlabel('n')
ax.set_ylabel('Overlapping Length (mm)')
ax.set_zlabel('Angle (degrees)')
ax.set_title('Angle as a function of n and Overlapping Length')

# Add colorbar
fig.colorbar(surf, label="Angle (degrees)")
plt.show()

