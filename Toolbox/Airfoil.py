from Airfoils import airfoilcoords
import numpy as np
from scipy.interpolate import interp1d
class Airfoil:
    def __init__(self, airfoilname, thickness, chordlength, Cl, moments, shearforces, sparlocations, topmembers, botmembers):
        """
        We make the assumption that the laminate does not change within a member!
        For now we'll assume the whole structure is sandwich panel?

        :param riblocations: locations from LE to TE of spars
        :param laminates: Members from rib to rib
        :param airfoilname:
        :param thickness:
        :param chordlength:
        :param Cl:
        :param moments:
        - list of moments -> there are 2, around both the x and y axis
        :param shearforces:
        - list of shear forces on the profile -> there are again 2, in the direction of both axes
        """
        self.airfoilname = airfoilname
        self.thickness = thickness
        self.chordlength = chordlength
        self.Cl = Cl

        self.topmembers = topmembers
        self.botmembers = botmembers

        # make coordinates:
        self.topcoordinates_base = airfoilcoords[airfoilname]['top']
        self.bottomcoordinates_base = airfoilcoords[airfoilname]['bot']
        self.top_bot_coordinates()

        self.sparlocations = sparlocations

        self.FindMemberLocations()

    def top_bot_coordinates(self):
        self.topcoordinates = self.topcoordinates_base
        self.botcoordinates = self.bottomcoordinates_base
        return

    def FindMemberLocations(self):
        # assign the width of the members based on the locations of the ribs
        for i, member in enumerate(self.botmembers):
            if i == 0:
                member.startcoord = [0, 0]
                member.endcoord = [self.sparlocations[i], self.Bot_height_at(self.sparlocations[i])]
            elif i == len(self.botmembers) - 1:
                # last member:
                member.startcoord = [self.sparlocations[i-1], self.Bot_height_at(self.sparlocations[i-1])]
                member.endcoord = [self.chordlength, 0]
            else:
                member.startcoord = [self.sparlocations[i-1], self.Bot_height_at(self.sparlocations[i-1])]
                member.endcoord = [self.sparlocations[i], self.Bot_height_at(self.sparlocations[i])]
            member.Calculate_b()

        for i, member in enumerate(self.topmembers):
            if i == 0:
                member.startcoord = [0, 0]
                member.endcoord = [self.sparlocations[i], self.Top_height_at(self.sparlocations[i])]
            elif i == len(self.botmembers) - 1:
                # last member:
                member.startcoord = [self.sparlocations[i-1], self.Top_height_at(self.sparlocations[i-1])]
                member.endcoord = [self.chordlength, 0]
            else:
                member.startcoord = [self.sparlocations[i-1], self.Top_height_at(self.sparlocations[i-1])]
                member.endcoord = [self.sparlocations[i], self.Top_height_at(self.sparlocations[i])]
            member.Calculate_b()
        return

    def Top_height_at(self, x):
        xlist = [point[0] for point in self.topcoordinates]
        ylist = [point[1] for point in self.topcoordinates]
        top_interp = interp1d(xlist, ylist, kind='linear')
        return top_interp(x)

    def Bot_height_at(self, x):
        xlist = [point[0] for point in self.botcoordinates]
        ylist = [point[1] for point in self.botcoordinates]
        bot_interp = interp1d(xlist, ylist, kind='linear')
        return bot_interp(x)

    def Cross_sectional_analysis(self):
        return

    def segment_I(self, Ex, t, P1, P2):
        """"
        Calculates Ixy, Ixx ,Iyy around neutral axis
        """
        # first Ixx:
        x1, y1 = P1
        x2, y2 = P2

        b = (y1 - (x1*y2/x2))/(1-x1/x2)
        a = (y2-b)/x2
        Ixx = t*np.sqrt(1+a**2)*(a**2*(x2**3-x1**3)/3 + a*b*(x2**2-x1**2) + b**2*(x2-x1))

        # Ixy in a similar way:
        Ixy = t*np.sqrt(1+a**2)*((a*x2**4)/4 - (a*x1**4)/4 + (b*x2**3)/3 - (b*x1**3)/3)

        # We can simply switch x and y, and obtain Iyy using the same equations:
        y1, x1 = P1
        y2, x2 = P2

        b = (y1 - (x1 * y2 / x2)) / (1 - x1 / x2)
        a = (y2 - b) / x2
        Iyy = t * np.sqrt(1 + a ** 2)*(a ** 2 * (x2 ** 3 - x1 ** 3) / 3 + a * b * (x2 ** 2 - x1 ** 2) + b ** 2 * (x2 - x1))

        return Ixx, Iyy, Ixy

    def Neutralpoints(self):

        return xneutral, yneutral

    def Curvatures(self):
        I = np.array([[Iyy, -Ixy],
                     [-Ixy, Ixx]])
