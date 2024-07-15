import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from Toolbox.Airfoil import Airfoil
from Toolbox.Member import Member
from Data.Panels import PanelWingRoot

class Wing:
    def __init__(self, liftdistribution, coorddistribution, LElocations, thicknessdistribution, halfspan, ribcoordinates, sparcoordinates):
        """
        :param airfoilcoords: list
        - list containing 2 lists
        - top airfoil coordinates
        - bottom airfoil coordinates

        :param halfspan: float
        - half span of the wing (on one side of the fuselage/center point)

        :param liftdistribution: list
        - list of lift per mm span along the wing, at any number of points, normalized to the half span
        - Constant interval of the points.
        - counts from the root to the tip

        :param coorddistribution: list
        - Constant interval of the points.
        - counts from the root to the tip
        - list of coord of the airfoil, in mm

        :param thicknessdistribution: list
        - Constant interval of the points.
        - counts from the root to the tip
        - list of thickness of the airfoil in % of the coord

        :param ribcoordinates: list
        - list containing the y coordinates along the half span of the wing

        :param sparcoordinates: list
        - list of one list per spar:
        -> containing sets of x, y coordinates for the points of the rib
        """

        self.liftdistribution = liftdistribution  # list
        self.chorddistribution = coorddistribution  # list
        self.LElocations = LElocations
        self.thicknessdistribution = thicknessdistribution # list
        self.halfspan = halfspan  # float
        self.sparcoordinates = sparcoordinates
        self.ribcoordinates = ribcoordinates
        self.toppanels = []
        self.botpanels = []

    def lift_at(self, y):
        """
        Interpolates the lift per unit span at a specific location along the span.
        :param location: float
        :return: float
        """
        listofpoints = np.linspace(0, self.halfspan, len(self.liftdistribution))
        lift_interp = interp1d(listofpoints, self.liftdistribution, kind='linear')
        return lift_interp(y)

    def LE_at(self, y):
        """
        Interpolates the leading edge coordinate
        :param location: float
        :return: float
        """
        listofpoints = np.linspace(0, self.halfspan, len(self.LElocations))
        LE_interp = interp1d(listofpoints, self.LElocations, kind='linear')
        return LE_interp(y)

    def chord_at(self, y):
        """
        Interpolates the leading edge coordinate
        :param location: float
        :return: float
        """
        listofpoints = np.linspace(0, self.halfspan, len(self.chorddistribution))
        LE_interp = interp1d(listofpoints, self.chorddistribution, kind='linear')
        return LE_interp(y)

    def thickness_at(self, y):
        """
        Interpolates the leading edge coordinate
        :param location: float
        :return: float
        """
        listofpoints = np.linspace(0, self.halfspan, len(self.thicknessdistribution))
        thickness_interp = interp1d(listofpoints, self.thicknessdistribution, kind='linear')
        return thickness_interp(y)

    def internal_moment(self):
        # Initialize the internal moment array
        moment_distribution_internal = np.zeros_like(self.liftdistribution)
        y = np.linspace(0, self.halfspan, len(self.liftdistribution))

        # For each spanwise location, integrate the lift from y to half-span
        for i in range(len(liftdistribution)):
            # Lift distribution from yi to b/2
            lift_segment = self.liftdistribution[i:]
            y_segment = y[i:]

            # Compute the moment arm (distance from yi to each segment point)
            moment_arm = y_segment - y[i]

            # Compute the internal moment
            moment_distribution_internal[i] = np.trapz(lift_segment * moment_arm, y_segment)

        self.moment_distribution_internal = moment_distribution_internal
        return moment_distribution_internal

    def shear_force(self):
        # Initialize the shear force array
        shear_distribution = np.zeros_like(self.liftdistribution)
        y = np.linspace(0, self.halfspan, len(self.liftdistribution))

        # For each spanwise location, integrate the lift from y to half-span
        for i in range(len(self.liftdistribution)):
            # Lift distribution from yi to b/2
            lift_segment = self.liftdistribution[i:]
            y_segment = y[i:]

            # Compute the shear force
            shear_distribution[i] = np.trapz(lift_segment, y_segment)

        self.shear_distribution = shear_distribution
        return shear_distribution

    def shear_at(self, location):
        """
        Interpolates the shear force at a specific location along the span.
        :param location: float
        :return: float
        """
        y = np.linspace(0, self.halfspan, len(self.shear_distribution))
        if not hasattr(self, 'shear_distribution'):
            self.shear_force()
        shear_interp = interp1d(y, self.shear_distribution, kind='linear')
        return shear_interp(location)

    def moment_at(self, location):
        """
        Interpolates the internal moment at a specific location along the span.
        :param location: float
        :return: float
        """
        y = np.linspace(0, self.halfspan, len(self.shear_distribution))
        if not hasattr(self, 'moment_distribution_internal'):
            self.internal_moment()
        moment_interp = interp1d(y, self.moment_distribution_internal, kind='linear')
        return moment_interp(location)

    def panels_at(self, objects_list, coordinate):
        for panel, end_coordinate in objects_list:
            if coordinate <= end_coordinate:
                return panel
        return None  # or raise an exception if you prefer

    def topmembers_at(self, y):
        # for a given location y, there must be members defined ->
        spars = self.spar_positions_at(y)
        panel = self.panels_at(self.toppanels, y)

        topmembers = [Member(panel) for _ in range(len(spars)+1)]
        return topmembers

    def botmembers_at(self, y):
        spars = self.spar_positions_at(y)
        panel = self.panels_at(self.botpanels, y)

        botmembers = [Member(panel) for _ in range(len(spars)+1)]
        return botmembers

    def curvatures_at(self, y):
        """
        Does analysis of a cross section of the wing at location y along the half span
        :return:
        """
        moment = self.moment_at(y)
        spars = self.spar_positions_at(y)-self.LE_at(y) # list of spar locations? or of objects? these indicate the locations of the spars for
                    # cross sectional analysis
        thickness = self.thickness_at(y)
        chordlength = self.chord_at(y)
        topmembers = self.topmembers_at(y)
        botmembers = self.botmembers_at(y)
        airfoil = Airfoil('e395-il', thickness, chordlength, spars, topmembers, botmembers)
        airfoil.Neutralpoints()
        airfoil.CalculateEI()
        kx, ky = airfoil.curvatures(moment, 0)
        return kx, ky

    def calculate_deflection(self, num_points=100):
        # Sample points along the half span
        y_points = np.linspace(0, self.halfspan, num=num_points)

        # Calculate curvatures at the sampled points
        curvatures = np.array([self.curvatures_at(y) for y in y_points])
        print(len(y_points), len(curvatures))

        # Integrate curvature to get slope
        slope = self.cumulative_trapezoid(curvatures, y_points, initial=0)

        # Integrate slope to get deflection
        deflection = self.cumulative_trapezoid(slope, y_points, initial=0)

        return y_points, deflection

    def cumulative_trapezoid(self, y, x, initial=0):
        """
        Perform cumulative trapezoidal integration by hand with an initial value.

        Parameters:
        y (array): The values to integrate.
        x (array): The x-coordinates corresponding to y.
        initial (float): The initial value to start the integration from.

        Returns:
        array: The cumulative integral of y with respect to x.
        """
        integral = np.zeros_like(y)
        integral[0] = initial
        for i in range(1, len(y)):
            integral[i] = integral[i - 1] + (y[i] + y[i - 1]) / 2 * (x[i] - x[i - 1])
        return integral

    def plot_deflection(self, num_points=100):
        y_points, deflection = self.calculate_deflection(num_points)

        # Plot the deflection
        plt.figure(figsize=(10, 6))
        plt.plot(y_points, deflection, label='Deflection')
        plt.xlabel('Half-span (y)')
        plt.ylabel('Deflection (v)')
        plt.title('Deflection of the Wing')
        plt.legend()
        plt.grid(True)
        plt.show()

    def spar_positions_at(self, location):
        """
        Returns the x coordinates of the spars at a given y location along the span.
        :param location: float
        :return: list of floats
        """
        spar_x_coords = []
        for spar in self.sparcoordinates:
            spar_x = []
            spar_y = []
            for coord in spar:
                spar_x.append(coord[0])
                spar_y.append(coord[1])
            spar_interp = interp1d(spar_y, spar_x, kind='linear', fill_value="extrapolate")
            spar_x_coords.append(spar_interp(location))
        return spar_x_coords

    def plot_wing_planform(self):
        """
        Plots the wing planform, including leading and trailing edges, spars, and ribs.
        :return: None
        """
        trailing_edge = [LE + chord for LE, chord in zip(self.LElocations, self.chorddistribution)]
        plt.figure(figsize=(12, 6))
        plt.plot(np.linspace(0, self.halfspan, len(self.LElocations)), self.LElocations, label='Leading Edge')
        plt.plot(np.linspace(0, self.halfspan, len(trailing_edge)), trailing_edge, label='Trailing Edge')
        plt.fill_between(np.linspace(0, self.halfspan, len(self.LElocations)), self.LElocations, trailing_edge, color='lightblue', alpha=0.5)

        # Plot the spars
        for spar in self.sparcoordinates:
            spar_x = [coord[0] for coord in spar]
            spar_y = [coord[1] for coord in spar]
            plt.plot(spar_y, spar_x, label='Spar', color='red')

        # Plot the ribs
        for rib in self.ribcoordinates:
            LE_at_rib = np.interp(rib, np.linspace(0, self.halfspan, len(self.LElocations)), self.LElocations)
            TE_at_rib = np.interp(rib, np.linspace(0, self.halfspan, len(trailing_edge)), trailing_edge)
            plt.plot([rib, rib], [LE_at_rib, TE_at_rib], color='green', linestyle='--', label='Rib' if rib == self.ribcoordinates[0] else "")

        plt.xlabel('Spanwise location (y)')
        plt.ylabel('Chordwise location (x)')
        plt.title('Wing Planform')
        plt.legend()
        plt.grid(True)
        plt.axis('equal')
        plt.show()

    def Failureanalysis(self):
        """"
        Does failure analysis of the wing
        :return:
        """
        return

def generate_chord_and_leading_edge(n, halfspan, coord_at_root):
    y = np.linspace(0, halfspan, n)
    chord_lengths = coord_at_root * np.sqrt(1 - (y / halfspan) ** 2)
    leading_edge_offsets = 0.25 * coord_at_root - 0.25 * chord_lengths
    leading_edge_locations = leading_edge_offsets + 0

    return list(chord_lengths), list(leading_edge_locations)

def generate_lift(n, halfspan, lift):
    y = np.linspace(0, halfspan, n)
    halflift = lift/2
    lr = (np.pi/4)*halflift/halfspan
    lifts = lr * np.sqrt(1 - (y / halfspan) ** 2)
    return lifts

# Example usage:
n = 1000  # Number of points
halfspan = 1500  # Half span of the wing
coord_at_root = 300  # Chord length at the root

chord_lengths, leading_edge_locations = generate_chord_and_leading_edge(n, halfspan, coord_at_root)

# Example usage:
liftdistribution = generate_lift(10, 1500, 981)
thicknessdistribution = [10, 12, 14, 16, 18]  # Example thickness distribution
halfspan = 1500  # Half span of the wing
sparcoordinates = [[[300/4, 0], [15 + 250/4, 1000], [75, 1500]], [[200, 0], [185, 1000], [75, 1500]]]  # Example spar coordinates
ribcoordinates = [0, 600, 1200, 1400]



wing = Wing(liftdistribution, chord_lengths, leading_edge_locations, thicknessdistribution, halfspan, ribcoordinates, sparcoordinates)

wing.toppanels = [[PanelWingRoot, 1600]]
wing.botpanels = [[PanelWingRoot, 1600]]
moment_distribution = wing.internal_moment()
shear_distribution = wing.shear_force()

location = 1200  # Example location along the half span
shear_at_location = wing.shear_at(location)
moment_at_location = wing.moment_at(location)
spar_positions = wing.spar_positions_at(location)
LE_at_location = wing.LE_at(location)

print("Shear Force at location", location, ":", shear_at_location)
print("Moment at location", location, ":", moment_at_location)
print("Spar positions at location", location, ":", spar_positions)
print("Leading Edge at location", location, ":", LE_at_location)
print('airfoil curvature:', location, ":", wing.curvatures_at(location))

# Plot the wing planform
wing.plot_wing_planform()
wing.plot_deflection(100)