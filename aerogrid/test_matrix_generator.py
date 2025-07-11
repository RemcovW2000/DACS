import csv
import numpy as np

n_range = np.linspace(0.1, 6, 5)
overlapping_range = np.linspace(40, 200, 5)
radius_range = np.linspace(8, 20, 7)

# Create the header row
header = ["Name", "P1", "P2", "P6", "P13","P14","P15","P16","P7", "P8", "P9", "P10"]

half_span = 1.5 # meters
weight = 25 # kg
g = 9.81 # m/s^2

def compute_loadcase(n):
    force = (weight * g /2) * n
    arm =  4 * half_span / (3 * np.pi)
    print(arm)
    arm_in_sim = 0.4 # m
    moment = (arm - arm_in_sim) * force
    return moment, -force

print(compute_loadcase(4))

# Generate the rows
rows = []
index = 0
for n in n_range:
    for o in overlapping_range:
        for r in radius_range:
            p2, p1 = compute_loadcase(n)
            p6 = 150 - o
            p13 = r - 1
            p14 = r
            p15 = r
            p16 = r + 1
            # Create a new row
            row = [f"DP {index}", p1, p2, p6, p13, p14, p15, p16, "", "", "", ""]  # P7 and P8 are empty
            rows.append(row)
            index += 1

# Write to a CSV file
with open("generated_data_01.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)  # Write the header
    writer.writerows(rows)   # Write all the rows

print("CSV file 'generated_data_01.csv' has been created.")