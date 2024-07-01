import numpy as np
from scipy.interpolate import interp1d


class Wing:
    def __init__(self, liftdistribution, coorddistribution, thicknessdistribution, halfspan, airfoilcoords, ribcoordinates, sparcoordinates):
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
        - list containing the y coordinates along the span of the locations of the ribs

        :param sparcoordinates: list
        - list containing the x coordinates along the coord of the wing
        """

        self.liftdistribution = np.array(liftdistribution)  # list
        self.coorddistribution = np.array(coorddistribution)  # list
        self.thicknessdistribution = np.array(thicknessdistribution)  # list
        self.halfspan = halfspan  # float
        self.n_points = len(self.liftdistribution)
        self.y = np.linspace(0, halfspan, self.n_points)  # spanwise locations
        self.dy = self.y[1] - self.y[0]  # interval between points (assumed constant)
        self.airfoilcoords = airfoilcoords

    def internal_moment(self):
        # Initialize the internal moment array
        moment_distribution_internal = np.zeros_like(self.liftdistribution)

        # For each spanwise location, integrate the lift from y to half-span
        for i in range(self.n_points):
            # Lift distribution from yi to b/2
            lift_segment = self.liftdistribution[i:]
            y_segment = self.y[i:]

            # Compute the moment arm (distance from yi to each segment point)
            moment_arm = y_segment - self.y[i]

            # Compute the internal moment
            moment_distribution_internal[i] = np.trapz(lift_segment * moment_arm, y_segment)

        self.moment_distribution_internal = moment_distribution_internal
        return moment_distribution_internal

    def shear_force(self):
        # Initialize the shear force array
        shear_distribution = np.zeros_like(self.liftdistribution)

        # For each spanwise location, integrate the lift from y to half-span
        for i in range(self.n_points):
            # Lift distribution from yi to b/2
            lift_segment = self.liftdistribution[i:]
            y_segment = self.y[i:]

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
        if not hasattr(self, 'shear_distribution'):
            self.shear_force()
        shear_interp = interp1d(self.y, self.shear_distribution, kind='linear')
        return shear_interp(location)

    def moment_at(self, location):
        """
        Interpolates the internal moment at a specific location along the span.
        :param location: float
        :return: float
        """
        if not hasattr(self, 'moment_distribution_internal'):
            self.internal_moment()
        moment_interp = interp1d(self.y, self.moment_distribution_internal, kind='linear')
        return moment_interp(location)

    def Cross_sectional_analysis(self, y):
        """
        Does analysis of a cross section of the wing at location y along the half span
        :return:
        """
        shear = self.shear_at(y)
        moment = self.moment_at(y)
        spars = [] # list of spar locations? or of objects? these indicate the locations of the spars for
                    # cross sectional analysis

        return

    def Generatepanels(self):
        """"
        Generates panels to analyze
        :return:
        """
        return

    def Failureanalysis(self):
        """"
        Does failure analysis of the wing
        :return:
        """
        return


# Example usage:
liftdistribution = [1, 2, 3, 4, 5]  # Example lift distribution
coorddistribution = [100, 200, 300, 400, 500]  # Example chord distribution
thicknessdistribution = [10, 12, 14, 16, 18]  # Example thickness distribution
halfspan = 5  # Half span of the wing

wing = Wing(liftdistribution, coorddistribution, thicknessdistribution, halfspan)
moment_distribution = wing.internal_moment()
shear_distribution = wing.shear_force()

location = 2.5  # Example location along the half span
shear_at_location = wing.shear_at(location)
moment_at_location = wing.moment_at(location)

print("Shear Force at location", location, ":", shear_at_location)
print("Moment at location", location, ":", moment_at_location)
