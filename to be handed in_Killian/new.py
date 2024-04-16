import numpy as np

list = np.array([1, 2])
list1 = np.array([1, 2])

if list.size > 0 and list1.size > 0:
    print('Both arrays are non-empty')
else:
    print('At least one array is empty')
