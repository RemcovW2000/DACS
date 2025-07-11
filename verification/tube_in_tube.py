import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load the CSV file (adjust path as needed)
file_path = "parameter_force_moment_result.csv"

# Load the file while skipping comment lines
df = pd.read_csv(file_path, comment="#", skip_blank_lines=True)

# Define the tube diameter
diameter = 0.030  # meters

# Compute the deformation angle (in degrees)
df["Angle (deg)"] = np.degrees(np.arctan((df["P7"] - df["P8"]) / diameter))
force_arm = 0.1 #meters
df["moment at joint"] = -df["P1"] * force_arm + df["P2"]
# Create a 3D scatter plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
print(df)
sc = ax.scatter(df["P1"], df["moment at joint"], df["Angle (deg)"], c=df["Angle (deg)"], cmap="viridis", edgecolors="k")

# Labels and title
ax.set_xlabel("Force X Component (N)")
ax.set_ylabel("Moment Y Component (N·m)")
ax.set_zlabel("Deformation Angle (deg)")
ax.set_title("Deformation Angle as Function of Force and Moment")

# Add color bar
cbar = plt.colorbar(sc, ax=ax, shrink=0.5, aspect=10)
cbar.set_label("Deformation Angle (deg)")

plt.show()