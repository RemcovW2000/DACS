import numpy as np

# Create a sample 6x4 numpy array
array = np.array([[1, 2, 3, 4],
                  [5, 6, 7, 8],
                  [9, 10, 11, 12],
                  [13, 14, 15, 16],
                  [17, 18, 19, 20],
                  [21, 22, 23, 24]])

# Extract elements from the 2nd and 3rd rows for each column
tuples = tuple(zip(array[1], array[2]))
print(tuples)
print(tuples[2])
