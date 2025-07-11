import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm

# Function to load data from file
def load_data(filename):
    data = np.loadtxt(filename, skiprows=1)  # Skip the header
    x = data[:, 1]  # X Coordinate
    y = data[:, 2]  # Y Coordinate
    shear_stress = data[:, 4]  # Shear Stress
    return x, y, shear_stress

# Load data from files
rear_x, rear_y, rear_xz = load_data("TS_xz.txt")
rear_x, rear_y, rear_yz = load_data("TS_yz.txt")

front_x, front_y, front_xz = load_data("TS_xz.txt")
front_x, front_y, front_yz = load_data("TS_yz.txt")

airfoil_x, airfoil_y, airfoil_xz = load_data("TS_xz.txt")
airfoil_x, airfoil_y, airfoil_yz = load_data("TS_yz.txt")

class Segment:
    def __init__(self, p1, p2, direction, magnitude, force_direction):
        self.p1 = p1
        self.p2 = p2
        self.direction = direction
        self.magnitude = magnitude
        self.force_direction = force_direction


def compute_segments(x, y, xz, yz):
    segments = []
    for i in range(len(x) - 1):
        p1 = (x[i], y[i])
        p2 = (x[i + 1], y[i + 1])
        vector = np.array([x[i + 1] - x[i], y[i + 1] - y[i]])
        vector /= np.linalg.norm(vector)
        force_vector = np.array([xz[i], yz[i]])
        force_vector /= np.linalg.norm(force_vector)
        magnitude = np.linalg.norm([xz[i], yz[i]])
        if np.dot(vector, force_vector) < 0:
            magnitude *= -1
        segments.append(Segment(p1, p2, vector, magnitude, force_vector))

    # Last point's vector is opposite of the previous
    if len(x) > 1:
        last_p1 = (x[-2], y[-2])
        last_p2 = (x[-1], y[-1])
        last_vector = -segments[-1].direction
        last_force_vector = np.array([xz[-1], yz[-1]])
        last_force_vector /= np.linalg.norm(last_force_vector)
        last_magnitude = np.linalg.norm([xz[-1], yz[-1]])
        if np.dot(last_vector, last_force_vector) < 0:
            last_magnitude *= -1
        segments.append(Segment(last_p1, last_p2, last_vector, last_magnitude, last_force_vector))

    return segments


def plot_segment_magnitude(segments, title):
    values = [seg.magnitude for seg in segments]
    min_value, max_value = min(values), max(values)
    max_abs_value = max(abs(min_value), abs(max_value))

    if max_value <= 0:
        cmap = plt.get_cmap('Blues')
        norm = Normalize(vmin=0, vmax=abs(min_value))
        normalized_values = [abs(val) for val in values]
    elif min_value >= 0:
        cmap = plt.get_cmap('Reds')
        norm = Normalize(vmin=0, vmax=max_value)
        normalized_values = values
    else:
        cmap = plt.get_cmap('coolwarm')
        norm = TwoSlopeNorm(vmin=-max_abs_value, vcenter=0, vmax=max_abs_value)
        normalized_values = values

    fig, ax = plt.subplots()
    for seg, value in zip(segments, normalized_values):
        x1, y1 = seg.p1
        x2, y2 = seg.p2
        ax.plot([x1, x2], [y1, y2], color=cmap(norm(value)), linewidth=2)

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_title(title)
    ax.grid()
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label='Magnitude')
    plt.show()


# Load data
rear_x, rear_y, rear_xz = load_data("TS_xz.txt")
rear_x, rear_y, rear_yz = load_data("TS_yz.txt")
rear_segments = compute_segments(rear_x, rear_y, rear_xz, rear_yz)
plot_segment_magnitude(rear_segments, "Rear Segment Magnitudes")

front_x, front_y, front_xz = load_data("FS_xz.txt")
front_x, front_y, front_yz = load_data("FS_yz.txt")
front_segments = compute_segments(front_x, front_y, front_xz, front_yz)
plot_segment_magnitude(front_segments, "Front Segment Magnitudes")

airfoil_x, airfoil_y, airfoil_xz = load_data("airfoil_xz.txt")
airfoil_x, airfoil_y, airfoil_yz = load_data("airfoil_yz.txt")
print(airfoil_x)
airfoil_segments = compute_segments(airfoil_x, airfoil_y, airfoil_xz, airfoil_yz)

all_segments = rear_segments + front_segments + airfoil_segments
plot_segment_magnitude(all_segments, "Airfoil Segment Magnitudes")
