import numpy as np
from MonteCarloSimulation import MonteCarloSimulation
# Assuming the mean properties and standard deviations for the materials are given as follows:
mean_properties = {
    'E1': 145300, 'E2': 8500, 'G12': 4580, 'v12': 0.31,
    'Xt': 1932, 'Xc': 1480, 'Yt': 108, 'Yc': 220, 'S': 132.8 ,
    'E11f': 500000, 'msf': 1.1, 'v21f': 0.1
}

std_dev_properties = {
    'E1': 3280, 'E2': 1280 , 'G12': 830, 'v12': 0.018,
    'Xt': 128.3, 'Xc': 100, 'Yt': 8.2, 'Yc': 17, 'S': 6.21,
    'E11f': 0, 'msf': 0, 'v21f': 0
}

# Configuration for a [0/90/45/-45]s laminate
laminas_configuration = [
(0.25, 0, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, 90, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, 45, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, -45, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, -45, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, 45, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, 90, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, 0, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, 0, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, 90, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, 45, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, -45, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, -45, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, 45, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, 90, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
(0.25, 0, {'elastic': ['E1', 'E2', 'G12', 'v12'], 'failure': ['E11f', 'v21f', 'msf', 'Xt', 'Xc', 'Yt', 'Yc', 'S']}),
]

# The load we apply to the simulation:
load = np.array([[500], [0], [0], [0], [0], [0]])
# load = np.array([[1932*2], [0], [0], [0], [0], [0]])

# Number of simulations to run
num_simulations = 1000

# Instantiate and run the Monte Carlo simulation
mc_sim = MonteCarloSimulation(num_simulations, mean_properties, std_dev_properties, laminas_configuration)
mc_sim.run(load)

# Analyze the results
mc_sim.analyze_results()
