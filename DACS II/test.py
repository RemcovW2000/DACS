import numpy as np

# Number of rows
n = 10

# Create a linspace array for the first column
linspace_array = np.linspace(0, 1, n)

# Initialize an array of zeros with n rows and 3 columns
data_array = np.zeros((n, 3))

# Fill the first column with the linspace array
data_array[:, 0] = linspace_array

list = []
print(len(list))