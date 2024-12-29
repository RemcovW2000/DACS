import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from Toolbox.airfoil import Airfoil
from Toolbox.member import Member
from Data.Panels import Sandwiches, Laminates
from scipy import integrate
from structural_entity import StructuralEntity

class Wing(StructuralEntity):
    """
    Wing Class
    """
    def __init__(self, lift_distribution, chord_distribution, LE_locations, thickness_distribution, half_span, rib_coordinates, spar_coordinates):
        """
        :param airfoilcoords: list
        - list containing 2 lists
        - top airfoil coordinates
        - bottom airfoil coordinates

        :param half_span: float
        - half span of the wing (on one side of the fuselage/center point)

        :param lift_distribution: list
        - list of lift per mm span along the wing, at any number of points, normalized to the half span
        - Constant interval of the points.
        - counts from the root to the tip

        :param chord_distribution: list
        - Constant interval of the points.
        - counts from the root to the tip
        - list of coord of the airfoil, in mm

        :param thickness_distribution: list
        - Constant interval of the points.
        - counts from the root to the tip
        - list of thickness of the airfoil in % of the coord

        :param rib_coordinates: list
        - list containing the y coordinates along the half span of the wing

        :param spar_coordinates: list
        - list of one list per spar:
        -> containing sets of x, y coordinates for the points of the rib
        """
        super().__init__('wing')
        self.lift_distribution = lift_distribution            # List containing lift per unit span as function of partial span coordinate
        self.chord_distribution = chord_distribution          # List containing chord as function of partial span coordinate y
        self.LE_locations = LE_locations                      # List containing the location of the leading edge as a function of partial span coordinate y
        self.thickness_distribution = thickness_distribution  # list defining the thickness of the wing as a function of the partial span coordinate y, NOT CURRENLY USED!
        self.half_span = half_span                            # The half span of the wing
        self.spar_coordinates = spar_coordinates              # List containing coordinates that define the location of the spars
        self.rib_coordinates = rib_coordinates                # Y coordinates of ribs
        self.toppanels = []     # Panels used for the top skin, as function of y coordinate
        self.botpanels = []     # Panels used for the bot skin, as function of y coordinate
        self.sparpanels = []    # Panels used for the spars, as function of y coordinate

        self.trcoordinates = [] # List containing the coordinates that define the location of the top reinforcement
        self.brcoordinates = [] # List containing the coordinates that define the location of the bottom reinforcement
        self.trwidth = None     # Width of top reinforcement
        self.brwidth = None     # Width of bottom reinforcement

        self.weight = None

        self.trpanels = []      # Panel used for top reinforcement
        self.brpanels = []      # Panel used for bottom reinforcement

        self.tip_buffer = 50   # Distance from the tip at which we no longer analyze the wing due to (close to) 0 thickness and load
        self.airfoils = []      # list of airfoil objects along half span of wing
        self.airfoil_FIs = []
        self.nr_airfoils = 10   # number of airfoils places in the cross-sectional analysis

    @property
    def child_objects(self):
        return self.airfoils

    def lift_at(self, y):
        """
        Interpolates the lift per unit span at a specific location along the span.
        :param location: float
        :return: float
        """
        listofpoints = np.linspace(0, self.half_span, len(self.lift_distribution))
        lift_interp = interp1d(listofpoints, self.lift_distribution, kind='linear')
        return lift_interp(y)

    def LE_at(self, y):
        """
        Interpolates the leading edge coordinate
        :param location: float
        :return: float
        """
        listofpoints = np.linspace(0, self.half_span, len(self.LE_locations))
        LE_interp = interp1d(listofpoints, self.LE_locations, kind='linear')
        return LE_interp(y)

    def chord_at(self, y):
        """
        Interpolates the chord lenght of the airoil
        :param location: float
        :return: float
        """
        listofpoints = np.linspace(0, self.half_span, len(self.chord_distribution))
        LE_interp = interp1d(listofpoints, self.chord_distribution, kind='linear')
        return LE_interp(y)

    def thickness_at(self, y):
        """
        NOT CURRENLY USED
        Interpolates the thickness of the airfoil
        :param location: float
        :return: float
        """
        listofpoints = np.linspace(0, self.half_span, len(self.thickness_distribution))
        thickness_interp = interp1d(listofpoints, self.thickness_distribution, kind='linear')
        return thickness_interp(y)

    def internal_moment(self):
        """
        Calculates the internal moment distribution in the wing as a function of the partial span coordinate y, due to
        external load.
        :return: list
        """
        # Initialize the internal moment array
        moment_distribution_internal = np.zeros_like(self.lift_distribution)
        y = np.linspace(0, self.half_span, len(self.lift_distribution))

        # For each spanwise location, integrate the lift from y to half-span
        for i in range(len(self.lift_distribution)):
            # Lift distribution from yi to b/2
            lift_segment = self.lift_distribution[i:]
            y_segment = y[i:]

            # Compute the moment arm (distance from yi to each segment point)
            moment_arm = y_segment - y[i]

            # Compute the internal moment
            moment_distribution_internal[i] = np.trapz(lift_segment * moment_arm, y_segment)

        self.moment_distribution_internal = moment_distribution_internal
        return moment_distribution_internal

    def shear_force(self):
        """
        Calculates the internal shear force distribution in the wing as a function of the partial span coordinate y, due
        to external load.
        :return: list
        """
        # Initialize the shear force array
        shear_distribution = np.zeros_like(self.lift_distribution)
        y = np.linspace(0, self.half_span, len(self.lift_distribution))

        # For each spanwise location, integrate the lift from y to half-span
        for i in range(len(self.lift_distribution)):
            # Lift distribution from yi to b/2
            lift_segment = self.lift_distribution[i:]
            y_segment = y[i:]

            # Compute the shear force
            shear_distribution[i] = np.trapz(lift_segment, y_segment)

        self.shear_distribution = shear_distribution
        return shear_distribution

    def shear_location_at(self, y):
        """
        Should find the point of application of the shear force in the cross section. Now just returns the x position of
        quarter chord, decent assumption.
        :return: float
        """
        return self.chord_at(0)/4, 0

    def shear_at(self, location):
        """
        Interpolates the shear force at a specific location along the span.
        :param location: float
        :return: float
        """
        y = np.linspace(0, self.half_span, len(self.shear_distribution))
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
        y = np.linspace(0, self.half_span, len(self.shear_distribution))
        if not hasattr(self, 'moment_distribution_internal'):
            self.internal_moment()
        moment_interp = interp1d(y, self.moment_distribution_internal, kind='linear')
        return moment_interp(location)

    def reinforcement_panel_at(self, location, side):
        """
        Finds the panel used for the reinforcement as a function of the y coordinate
        :param location: float
        :param side: string
        :return: panel object
        """
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

    def tr_start_end_at(self, location):
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

    def br_start_end_at(self, location):
        """
        Returns the x coordinate of the bot reinforcement at a given y location along the span.
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

    def generate_airfoil_objects(self):
        """
        This function generates the airfoil objects and records their y coordinates.

        It assigns the members and reinforcement correctly
        :return:
        """
        analysislocations = np.linspace(0, self.half_span - self.tip_buffer, self.nr_airfoils)
        # TODO: add twist in the frame of reference

        reinforcement_end = self.trpanels[-1][1]
        airfoils = []
        for y in analysislocations:
            top_members = self.top_members_at(y)
            bot_members = self.bot_members_at(y)
            spar_members = self.spar_members_at(y)
            reinforcementpaneltop = self.reinforcement_panel_at(y, 'top')
            reinforcementpanelbot = self.reinforcement_panel_at(y, 'bot')
            spar_locations = [self.spar_positions_at(y)[i] - self.LE_at(y) for i in range(len(self.spar_positions_at(y)))]

            if y <= reinforcement_end:
                trstart, trend = self.tr_start_end_at(y) - self.LE_at(y)
                brstart, brend = self.br_start_end_at(y) - self.LE_at(y)
            else:
                trstart, trend = None, None
                brstart, brend = None, None

            chord = self.chord_at(y)
            type = 'NACA2410'
            # make the option:
            airfoil = Airfoil(type, 1, chord, spar_locations, top_members, bot_members, spar_members, reinforcementpaneltop, trstart,
                              trend, reinforcementpanelbot, brstart, brend)
            airfoil.xshear = self.shear_location_at(y)[0] - self.LE_at(y)   # x coordinate of point of application in airfoil FOR
            airfoil.yshear = self.shear_location_at(y)[1]                   # y coordinate of point of application in airfoil FOR, without twist is equal to 0
            airfoils.append(airfoil)
            airfoil.y = y

        self.airfoils = airfoils
        return

    def find_curvatures(self):
        for airfoil in self.airfoils:
            airfoil.curvatures(airfoil.y)

    def solve_stresses_CSA(self):
        for airfoil in self.airfoils:
            y = airfoil.y
            mx = self.moment_at(y)
            my = 0
            sx = 0
            sy = self.shear_at(y)
            moments = [mx, my]
            shears = [sx, sy]
            center = [self.chord_at(y)/4, 0] # place it at quarter chord
            airfoil.solve_stresses_CSA(moments, shears, center)
        return

    def calculate_weight(self):
        airfoil_weights = [airfoil.calculate_weight_per_b() for airfoil in self.airfoils]
        airfoil_ys = [airfoil.y for airfoil in self.airfoils]
        self.weight = integrate.simps(airfoil_weights, airfoil_ys)
        return self.weight

    def panels_at(self, objects_list, coordinate):
        for panel, end_coordinate in objects_list:
            if coordinate <= end_coordinate:
                paneltoreturn = panel
            elif coordinate > end_coordinate:
                paneltoreturn = None
        return paneltoreturn

    def top_members_at(self, y):
        # for a given location y, there must be members defined ->
        spars = self.spar_positions_at(y)
        panel = copy.deepcopy(self.panels_at(self.toppanels, y))

        top_members = [Member(panel) for _ in range(len(spars)+1)]
        return top_members

    def bot_members_at(self, y):
        spars = self.spar_positions_at(y)
        panel = copy.deepcopy(self.panels_at(self.botpanels, y))

        bot_members = [Member(panel) for _ in range(len(spars)+1)]
        return bot_members

    def spar_members_at(self, y):
        spars = self.spar_positions_at(y)
        panel = copy.deepcopy(self.sparpanels)
        spar_members = [Member(panel) for _ in range(len(spars))]
        return spar_members

    def curvatures_at(self, y):
        """
        Does analysis of a cross section of the wing at location y along the half span
        :return:
        """
        moment = self.moment_at(y)
        spars = self.spar_positions_at(y)-self.LE_at(y) # list of spar locations? or of objects? these indicate the locations of the spars for
                    # cross sectional analysis
        thickness = self.thickness_at(y)
        chord_length = self.chord_at(y)
        top_members = self.top_members_at(y)
        bot_members = self.bot_members_at(y)
        spar_members = self.spar_members_at(y)
        airfoil = Airfoil('NACA2410', thickness, chord_length, spars, top_members, bot_members, spar_members, )
        airfoil.neutral_points()
        airfoil.calculate_EI()
        kx, ky = airfoil.curvatures(moment, 0)
        return kx, ky

    def calculate_deflection(self, num_points=10):
        # Sample points along the half span
        y_points = np.linspace(0, self.half_span, num=num_points)

        # Calculate curvatures at the sampled points
        curvatures = np.array([self.curvatures_at(y) for y in y_points])

        # Integrate curvature to get slope
        slope = self.cumulative_trapezoid(curvatures, y_points, initial=0)

        # Integrate slope to get deflection
        deflection = self.cumulative_trapezoid(slope, y_points, initial=0)

        return y_points, deflection

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
        for spar in self.spar_coordinates:
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
        trailing_edge = [LE + chord for LE, chord in zip(self.LE_locations, self.chord_distribution)]
        spanwise_locations = np.linspace(0, self.half_span, len(self.LE_locations))

        plt.figure(figsize=(9, 12))

        # Swap x and y for leading edge and trailing edge
        plt.plot(self.LE_locations, spanwise_locations, label='Leading Edge')
        plt.plot(trailing_edge, spanwise_locations, label='Trailing Edge')
        plt.fill_betweenx(spanwise_locations, self.LE_locations, trailing_edge, color='lightblue', alpha=0.5)

        # Plot the spars with x and y swapped
        for spar in self.spar_coordinates:
            spar_x = [coord[0] for coord in spar]
            spar_y = [coord[1] for coord in spar]
            plt.plot(spar_x, spar_y, label='Spar', color='red')

        # Plot the ribs with x and y swapped
        for rib in self.rib_coordinates:
            LE_at_rib = np.interp(rib, spanwise_locations, self.LE_locations)
            TE_at_rib = np.interp(rib, spanwise_locations, trailing_edge)
            plt.plot([LE_at_rib, TE_at_rib], [rib, rib], color='green', linestyle='--',
                     label='Rib' if rib == self.rib_coordinates[0] else "")

        # Plot the reinforcement:
        if self.trcoordinates:
            tr_x_front, tr_x_back, tr_y = self.get_reinforcement_coordinates(self.trcoordinates, self.trpanels, self.trwidth)

            plt.plot(tr_x_front, tr_y, color='black', label='Top Reinforcement Front Edge')
            plt.plot(tr_x_back, tr_y, color='black', label='Top Reinforcement Back Edge')
            plt.plot([tr_x_front[-1], tr_x_back[-1]], [tr_y[-1], tr_y[-1]], color='black')

            plt.fill_betweenx(tr_y, tr_x_front, tr_x_back, color='black', alpha=0.5)

        # Plot the bottom reinforcements
        if self.brcoordinates:
            br_x_front, br_x_back, br_y = self.get_reinforcement_coordinates(self.brcoordinates, self.brpanels, self.brwidth)

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

    def get_reinforcement_coordinates(self, coordinates, panels, width):
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

    def failure_analysis(self):
        """"
        Executes the failure_analysis() method in each airfoil object in self.airfoils
        :return:
        """
        self.airfoil_FIs = [airfoil.failure_analysis() for airfoil in self.airfoils]

        self.finalize_failure_analysis(None)

        return max(value for key, value in self.failure_indicators.items() if isinstance(value, (int, float)))

    def plot_max_FI(self):
        # Plot the deflection
        plt.figure(figsize=(10, 6))
        plt.plot([airfoil.y for airfoil in self.airfoils], self.airfoil_FIs, label='FI_max')
        plt.xlabel('Half-span (y)')
        plt.ylabel('Max failure indicator')
        plt.title('Maximum failure indicator at each airfoil')
        plt.legend()
        plt.grid(True)
        return plt

    @staticmethod
    def cumulative_trapezoid(y, x, initial=0):
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

def generate_chord_and_leading_edge_elyptical(n, half_span, coord_at_root):
    y = np.linspace(0, half_span, n)
    chord_lengths = coord_at_root * np.sqrt(1 - (y / half_span) ** 2)
    leading_edge_offsets = 0.25 * coord_at_root - 0.25 * chord_lengths
    leading_edge_locations = leading_edge_offsets + 0

    return list(chord_lengths), list(leading_edge_locations)

def generate_chord_and_leading_edge_tapered(n, half_span, chord_at_root, taper_ratio):
    ys_for_interp = [0, half_span]
    tip_chord = taper_ratio * chord_at_root
    chords_for_interp = [chord_at_root, tip_chord]

    ys = np.linspace(0, half_span, n)

    chords = [np.interp(y_now, ys_for_interp, chords_for_interp) for y_now in ys]

    leading_edge_locations = [0.25*(chord_at_root - chord) for chord in chords]

    return list(chords), list(leading_edge_locations)
def generate_lift_elyptical(n, half_span, lift):
    y = np.linspace(0, half_span, n)
    halflift = lift/2
    lr = (np.pi/4)*halflift/half_span
    lifts = lr * np.sqrt(1 - (y / half_span) ** 2)
    return lifts