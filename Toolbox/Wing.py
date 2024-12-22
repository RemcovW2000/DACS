import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from Toolbox.Airfoil import Airfoil
from Toolbox.Member import Member
from Data.Panels import Sandwiches, Laminates

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
        self.sparpanels = []

        self.trcoordinates = []
        self.brcoordinates = []
        self.trwidth = None
        self.brwidth = None

        self.weight = 0.5       # Total weight of wing in kg

        # TODO: put in some sort of global top- and bottom reinforcement panel, could be in same format as normal panels
        self.trpanels = []
        self.brpanels = []
        # TODO: put in a global way to find the width of the reinforcement

        self.tip_buffer = 200   # Distance from the tip at which we no longer analyze the wing due to (close to) 0 thickness and load
        self.airfoils = []      # list of airfoil objects along half span of wing
        self.airfoilys = []     # Y locations of airfoil sections to be simulated
        self.airoil_maxFIs = [] # Max failure indicator in airfoil sections

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
        for i in range(len(self.liftdistribution)):
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

    def shear_location_at(self, y):
        return self.chord_at(0)/4, 0

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

    def reinforcementpanel_at(self, location, side):
        # for a given location y, members are defined
        if side == 'top':
            panels = self.trpanels
        elif side == 'bot':
            panels = self.brpanels

        panel = self.panels_at(panels, location)
        if panel:
            panel = copy.deepcopy(panel)
            rpanel = Member(panel)
            return rpanel
        else:
            return None

    def trstartend_at(self, location):
        """
        Returns the x coordinate of the top reinforcement at a given y location along the span.
        :param location: float
        :return: list of floats
        """
        tr_x = []
        tr_y = []
        reinforcement_end = self.trpanels[-1][1]
        if location > reinforcement_end:
            return None, None
        for coord in self.trcoordinates:
            tr_x.append(coord[0])
            tr_y.append(coord[1])
        tr_interp = interp1d(tr_y, tr_x, kind='linear', fill_value="extrapolate")
        return tr_interp(location) - self.trwidth/2, tr_interp(location) + self.trwidth/2

    def brstartend_at(self, location):
        """
        Returns the x coordinate of the top reinforcement at a given y location along the span.
        :param location: float
        :return: list of floats
        """
        br_x = []
        br_y = []
        reinforcement_end = self.trpanels[-1][1]
        if location > reinforcement_end:
            return None, None
        for coord in self.brcoordinates:
            br_x.append(coord[0])
            br_y.append(coord[1])
        tr_interp = interp1d(br_y, br_x, kind='linear', fill_value="extrapolate")
        return tr_interp(location) - self.brwidth/2, tr_interp(location) + self.brwidth/2

    def GenerateAirfoils(self):
        """
        This function generates the airfoil objects and records their coordinates.

        It assigns the members and reinforcement correctly
        :return:
        """
        analysislocations = np.linspace(0, self.halfspan - self.tip_buffer, 10)
        # TODO: add twist in the frame of reference

        airfoils = []
        for y in analysislocations:
            topmembers = self.topmembers_at(y)
            botmembers = self.botmembers_at(y)
            sparmembers = self.sparmembers_at(y)
            reinforcementpaneltop = self.reinforcementpanel_at(y, 'top')
            reinforcementpanelbot = self.reinforcementpanel_at(y, 'bot')
            sparlocations = [self.spar_positions_at(y)[i] - self.LE_at(y) for i in range(len(self.spar_positions_at(y)))]

            trstart, trend = self.trstartend_at(y)
            brstart, brend = self.brstartend_at(y)
            chord = self.chord_at(y)
            type = 'NACA2410'
            # make the option:
            airfoil = Airfoil(type, 1, chord, sparlocations, topmembers, botmembers, sparmembers, reinforcementpaneltop, trstart,
                              trend, reinforcementpanelbot, brstart, brend)
            airfoil.xshear = self.shear_location_at(y)[0] - self.LE_at(y)   # x coordinate of point of application in airfoil FOR
            airfoil.yshear = self.shear_location_at(y)[1]                   # y coordinate of point of application in airfoil FOR, without twist is equal to 0
            airfoils.append(airfoil)

        self.airfoils = airfoils
        self.airfoilys = analysislocations
        return

    def findcurvatures(self):
        for i, airfoil in enumerate(self.airfoils):
            y = self.airfoilys[i]
            airfoil.curvatures(y)
        return

    def SolveStresses_CSA(self):
        for i, airfoil in enumerate(self.airfoils):
            y = self.airfoilys[i]
            mx = self.moment_at(y)
            my = 0
            sx = 0
            sy = self.shear_at(y)
            moments = [mx, my]
            shears = [sx, sy]
            center = [self.chord_at(y)/4, 0] # place it at quarter chord
            airfoil.SolveStresses_CSA(moments, shears, center)
        return

    def panels_at(self, objects_list, coordinate):
        for panel, end_coordinate in objects_list:
            if coordinate <= end_coordinate:
                paneltoreturn = panel
            elif coordinate > end_coordinate:
                paneltoreturn = None
        return paneltoreturn

    def topmembers_at(self, y):
        # for a given location y, there must be members defined ->
        spars = self.spar_positions_at(y)
        panel = copy.deepcopy(self.panels_at(self.toppanels, y))

        topmembers = [Member(panel) for _ in range(len(spars)+1)]
        return topmembers

    def botmembers_at(self, y):
        spars = self.spar_positions_at(y)
        panel = copy.deepcopy(self.panels_at(self.botpanels, y))

        botmembers = [Member(panel) for _ in range(len(spars)+1)]
        return botmembers

    def sparmembers_at(self, y):
        spars = self.spar_positions_at(y)
        panel = copy.deepcopy(self.sparpanels)
        sparmembers = [Member(panel) for _ in range(len(spars))]
        return sparmembers

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
        sparmembers = self.sparmembers_at(y)
        airfoil = Airfoil('NACA2410', thickness, chordlength, spars, topmembers, botmembers, sparmembers, )
        airfoil.Neutralpoints()
        airfoil.CalculateEI()
        kx, ky = airfoil.curvatures(moment, 0)
        return kx, ky

    def calculate_deflection(self, num_points=10):
        # Sample points along the half span
        y_points = np.linspace(0, self.halfspan, num=num_points)

        # Calculate curvatures at the sampled points
        curvatures = np.array([self.curvatures_at(y) for y in y_points])

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
        Plots the wing planform, including leading and trailing edges, spars, ribs, and reinforcements, with x and y axes swapped.
        :return: None
        """
        trailing_edge = [LE + chord for LE, chord in zip(self.LElocations, self.chorddistribution)]
        spanwise_locations = np.linspace(0, self.halfspan, len(self.LElocations))

        plt.figure(figsize=(9, 12))

        # Swap x and y for leading edge and trailing edge
        plt.plot(self.LElocations, spanwise_locations, label='Leading Edge')
        plt.plot(trailing_edge, spanwise_locations, label='Trailing Edge')
        plt.fill_betweenx(spanwise_locations, self.LElocations, trailing_edge, color='lightblue', alpha=0.5)

        # Plot the spars with x and y swapped
        for spar in self.sparcoordinates:
            spar_x = [coord[0] for coord in spar]
            spar_y = [coord[1] for coord in spar]
            plt.plot(spar_x, spar_y, label='Spar', color='red')

        # Plot the ribs with x and y swapped
        for rib in self.ribcoordinates:
            LE_at_rib = np.interp(rib, spanwise_locations, self.LElocations)
            TE_at_rib = np.interp(rib, spanwise_locations, trailing_edge)
            plt.plot([LE_at_rib, TE_at_rib], [rib, rib], color='green', linestyle='--',
                     label='Rib' if rib == self.ribcoordinates[0] else "")

        # Plot the reinforcement:
        if self.trcoordinates:
            tr_x_front, tr_x_back, tr_y = self.GetReinforcementsCoords(self.trcoordinates, self.trpanels, self.trwidth)

            plt.plot(tr_x_front, tr_y, color='black', label='Top Reinforcement Front Edge')
            plt.plot(tr_x_back, tr_y, color='black', label='Top Reinforcement Back Edge')
            plt.plot([tr_x_front[-1], tr_x_back[-1]], [tr_y[-1], tr_y[-1]], color='black')

            plt.fill_betweenx(tr_y, tr_x_front, tr_x_back, color='black', alpha=0.5)

        # Plot the bottom reinforcements
        if self.brcoordinates:
            br_x_front, br_x_back, br_y = self.GetReinforcementsCoords(self.brcoordinates, self.brpanels, self.brwidth)

            plt.plot(br_x_front, br_y, color='blue', label='bot Reinforcement Front Edge')
            plt.plot(br_x_back, br_y, color='blue', label='bot Reinforcement Back Edge')
            plt.plot([br_x_front[-1], br_x_back[-1]], [br_y[-1], br_y[-1]], color='blue')

            plt.fill_betweenx(br_y, br_x_front, br_x_back, color='blue', alpha=0.5)

        # plotting settings
        plt.xlabel('Chordwise location (x)')
        plt.ylabel('Spanwise location (y)')
        plt.title('Wing Planform')
        plt.legend()
        plt.grid(True)
        plt.axis('equal')

        return plt

    def GetReinforcementsCoords(self, coordinates, panels, width):
        # find reinforcement coordinates:
        tr_x = [coord[0] for coord in coordinates]
        tr_y = [coord[1] for coord in coordinates]

        # find end coordinate:
        y_end = panels[-1][-1]
        # find x coordinate at y coordinate:
        f = interp1d(tr_y, tr_x, kind='linear')
        x_end = f(y_end)

        # make new list:
        tr_y = [y for y in tr_y if y < y_end]
        tr_y.append(y_end)

        tr_x = [tr_x[i] for i in range(len(tr_y) - 1)]
        tr_x.append(x_end)

        tr_x_front = [tr_x[i] + width / 2 for i in range(len(tr_x))]
        tr_x_back = [tr_x[i] - width / 2 for i in range(len(tr_x))]
        return tr_x_front, tr_x_back, tr_y

    def Failureanalysis(self):
        """"
        Does failure analysis of the wing
        :return:
        """
        airfoil_maxFIs = []
        for airfoil in self.airfoils:
            airfoil_maxFI = airfoil.FailureAnalysis_CSA()
            airfoil_maxFIs.append(airfoil_maxFI)

        self.airoil_maxFIs = airfoil_maxFIs
        return

    def plot_maxFI(self):
        # Plot the deflection
        plt.figure(figsize=(10, 6))
        plt.plot(self.airfoilys, self.airoil_maxFIs, label='FI_max')
        plt.xlabel('Half-span (y)')
        plt.ylabel('Max failure indicator')
        plt.title('Maximum failure indicator at each airfoil')
        plt.legend()
        plt.grid(True)
        return plt

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