import pandas as pd
import numpy as np

# Load the original CSV explicitly as a comma-separated file
file_path = "parameter_length_force_moment.csv"
df = pd.read_csv(file_path, comment="#", skip_blank_lines=True, sep=",")

# Define new ranges for P1, P2, P3
p1_range = np.linspace(0, -200, 5)
p2_range = np.linspace(0, 200, 5)
p3_range = np.linspace(-100, 100, 5)

# Generate new rows
num = 1
new_rows = []
for p1 in p1_range:
    for p2 in p2_range:
        for p3 in p3_range:
            new_rows.append({
                "Name": f"DP {num}",
                "P1": p1,
                "P2": p2,
                "P3": '',
                "P4": '',
                "P6": p3,
            })
            num += 1

# Convert new rows to DataFrame
new_df = pd.DataFrame(new_rows)

# Append new rows to the original dataframe
df_modified = pd.concat([df, new_df], ignore_index=True)

# Save the modified CSV using comma as separator
modified_file_path = "modified_parameter_force_moment_result.csv"
df_modified.to_csv(modified_file_path, index=False, sep=",")

print(f"Modified CSV saved as {modified_file_path}")