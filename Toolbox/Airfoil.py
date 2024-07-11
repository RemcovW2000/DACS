import copy

from Airfoils import airfoilcoords
import numpy as np
from scipy.interpolate import interp1d
from Toolbox.Lamina import Lamina
from Toolbox import MP
from Toolbox.Laminate import Laminate, LaminateBuilder
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp
from Toolbox.DamagedRegion import *
from tqdm import tqdm
from Toolbox.Member import Member


class Airfoil:
    def __init__(self, airfoilname, thickness, chordlength, sparlocations, topmembers, botmembers):
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
        self.shearforces = shearforces

        self.topmembers = topmembers
        self.botmembers = botmembers

        # make coordinates:
        self.topcoordinates_base = airfoilcoords[airfoilname]['top']
        self.bottomcoordinates_base = airfoilcoords[airfoilname]['bottom']
        self.top_bot_coordinates()

        self.sparlocations = sparlocations

        self.FindMemberLocations()

    def top_bot_coordinates(self):
        # transform the coordinates as neccesary:
        topcoordinates = self.topcoordinates_base
        botcoordinates = self.bottomcoordinates_base

        topcoordinates = [[element * self.chordlength for element in sublist] for sublist in topcoordinates]
        botcoordinates = [[element * self.chordlength for element in sublist] for sublist in botcoordinates]
        self.topcoordinates = topcoordinates
        self.botcoordinates = botcoordinates
        return self.topcoordinates, self.botcoordinates

    def FindMemberLocations(self):
        # assign the width of the members based on the locations of the ribs
        for i, member in enumerate(self.botmembers):
            if i == 0:
                member.startcoord = [0, 0]
                member.endcoord = [self.sparlocations[i], self.Bot_height_at(self.sparlocations[i])]
            elif i == len(self.botmembers) - 1:
                # last member:
                member.startcoord = [self.sparlocations[i - 1], self.Bot_height_at(self.sparlocations[i - 1])]
                member.endcoord = [self.chordlength, 0]
            else:
                member.startcoord = [self.sparlocations[i - 1], self.Bot_height_at(self.sparlocations[i - 1])]
                member.endcoord = [self.sparlocations[i], self.Bot_height_at(self.sparlocations[i])]
            member.Calculate_b()

        for i, member in enumerate(self.topmembers):
            if i == 0:
                member.startcoord = [0, 0]
                member.endcoord = [self.sparlocations[i], self.Top_height_at(self.sparlocations[i])]
            elif i == len(self.botmembers) - 1:
                # last member:
                member.startcoord = [self.sparlocations[i - 1], self.Top_height_at(self.sparlocations[i - 1])]
                member.endcoord = [self.chordlength, 0]
            else:
                member.startcoord = [self.sparlocations[i - 1], self.Top_height_at(self.sparlocations[i - 1])]
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

    def segment_I(self, t, P1, P2):
        """"
        Calculates Ixy, Ixx ,Iyy around neutral axis
        """
        # first Ixx:
        x1, y1 = P1
        x2, y2 = P2

        b = (y1 - (x1 * y2 / x2)) / (1 - x1 / x2)
        a = (y2 - b) / x2
        Ixx = abs(t * np.sqrt(1 + a ** 2) * (
                    a ** 2 * (x2 ** 3 - x1 ** 3) / 3 + a * b * (x2 ** 2 - x1 ** 2) + b ** 2 * (x2 - x1)))

        # Ixy in a similar way:
        Ixy = t*np.sqrt(1+a**2) * ((x2**3 - x1**3)*a/3 + (x2**2-x1**2)*b/2)

        # We can simply switch x and y, and obtain Iyy using the same equations:
        y1, x1 = P1
        y2, x2 = P2

        if x1 == x2 or x2 == 0:
            x2 +=1e-10
            b = (y1 - (x1 * y2 / x2)) / (1 - x1 / x2)
            a = (y2 - b) / x2
        else:
            b = (y1 - (x1 * y2 / x2)) / (1 - x1 / x2)
            a = (y2 - b) / x2
        Iyy = abs(t * np.sqrt(1 + a ** 2) * (
                    a ** 2 * (x2 ** 3 - x1 ** 3) / 3 + a * b * (x2 ** 2 - x1 ** 2) + b ** 2 * (x2 - x1)))

        return Ixx, Iyy, Ixy

    def segment_q(self, t, P1, P2, Ex, shearforces):
        Sx, Sy = shearforces
        Ixx, Iyy, Ixy = self.segment_I(t, P1, P2)
        term11 = -(Sx*Ixx - Sy*Ixy)/(Ixx*Iyy - Ixy**2)
        term12 = 0
        term21 = -(Sy*Iyy - Sx*Ixy)/(Ixx*Iyy - Ixy**2)
        term22 = 0
        return


    def Neutralpoints(self):
        # find the neutral points -> x and y coordinate of neutral axes around x and y axes:
        EAytot = 0
        EAxtot = 0
        EAtot = 0
        for member in self.topmembers:
            npoints = 100  # TODO: placeholder!
            # in the topmembers list, we know the start and end coordinates:
            xlist = np.linspace(member.startcoord[0], member.endcoord[0], npoints)
            for i, x in enumerate(xlist):
                if i == len(xlist) - 1:
                    break
                elif x == 0:
                    p1 = [0, 0]
                    p2 = [xlist[i + 1], self.Top_height_at(xlist[i + 1])]
                else:
                    # we have 1 point now:
                    p1 = [x, self.Top_height_at(x)]
                    p2 = [xlist[i + 1], self.Top_height_at(xlist[i + 1])]

                # given these 2 points we can find EA eq and distance:
                Dl = np.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)
                EA = Dl * member.panel.h * member.panel.Ex
                Dx = (p1[0] + p2[0]) / 2
                Dy = (p1[1] + p2[1]) / 2

                EAy = EA * Dy
                EAx = EA * Dx

                EAytot += EAy
                EAxtot += EAx
                EAtot += EA

        for member in self.botmembers:
            npoints = 100  # TODO: placeholder!
            # in the topmembers list, we know the start and end coordinates:
            xlist = np.linspace(member.startcoord[0], member.endcoord[0], npoints)
            for i, x in enumerate(xlist):
                if i == len(xlist) - 1:
                    break
                elif x == 0:
                    p1 = [0, 0]
                    p2 = [xlist[i + 1], self.Bot_height_at(xlist[i + 1])]
                else:
                    # we have 1 point now:
                    p1 = [x, self.Bot_height_at(x)]
                    p2 = [xlist[i + 1], self.Bot_height_at(xlist[i + 1])]

                # given these 2 points we can find EA eq and distance:
                length = np.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)
                EA = length * member.panel.h * member.panel.Ex
                Dx = (p1[0] + p2[0]) / 2
                Dy = (p1[1] + p2[1]) / 2

                EAy = EA * Dy
                EAx = EA * Dx

                EAytot += EAy
                EAxtot += EAx
                EAtot += EA

        ybar = EAytot / EAtot
        xbar = EAxtot / EAtot

        self.xbar = xbar
        self.ybar = ybar
        return xbar, ybar

    def plotairfoil(self):
        xlist = np.linspace(0.1, self.chordlength, 100)
        top = [self.Top_height_at(x) for x in xlist]
        bot = [self.Bot_height_at(x) for x in xlist]


        # Plotting
        plt.figure(figsize=(10, 6))
        plt.plot(xlist, top, label='Top Height')
        plt.plot(xlist, bot, label='Bottom Height')
        plt.axhline(y=self.ybar, linestyle='--')
        plt.axvline(x=self.xbar, linestyle='--')
        plt.axis('equal')
        plt.xlabel('X values')
        plt.ylabel('Heights')
        plt.title('Top and Bottom Heights vs. X values')
        plt.legend()
        plt.grid(True)
        plt.show()

    def CalculateEI(self):
        EIxx, EIyy, EIxy = 0, 0, 0
        for member in self.topmembers:
            npoints = 100  # TODO: placeholder!
            # in the topmembers list, we know the start and end coordinates:
            xlist = np.linspace(member.startcoord[0], member.endcoord[0], npoints)
            for i, x in enumerate(xlist):
                if i == len(xlist) - 1:
                    break
                elif x == 0:
                    p1 = [0, 0]
                    p2 = [xlist[i + 1], self.Bot_height_at(xlist[i + 1])]
                else:
                    # we have 1 point now:
                    p1 = [x, self.Top_height_at(x)]
                    p2 = [xlist[i + 1], self.Top_height_at(xlist[i + 1])]

                Ilist = self.segment_I(member.panel.h, p1, p2)
                EIxx += Ilist[0]*member.panel.Ex
                EIyy += Ilist[1]*member.panel.Ex
                EIxy += Ilist[2]*member.panel.Ex
        self.EIxx = EIxx
        self.EIyy = EIyy
        self.EIxy = EIxy
        return EIxx, EIyy, EIxy

    def curvatures(self, Mx, My):
        EIxx = self.EIxx
        EIyy = self.EIyy
        EIxy = self.EIxy
        kx = (Mx * EIyy - My * EIxy) / (EIxx * EIyy - EIxy ** 2)
        print('kx', kx)
        ky = (My * EIxx - Mx * EIxy) / (EIxx * EIyy - EIxy ** 2)
        return kx, ky

    def AxialStrain(self, kx, ky, x, y):
        e = -kx*(y-self.ybar) + ky*(x-self.xbar)
        return e



Laminate = LaminateBuilder([0], False, True, 1)
Member = Member(Laminate)
sparlocations = [80, 200]
topmembers = [copy.deepcopy(Member) for _ in range(len(sparlocations))]
botmembers = [copy.deepcopy(Member) for _ in range(len(sparlocations))]

type = 'e395-il'
Airfoil = Airfoil(type, 1, 300, None, sparlocations, topmembers, botmembers)
Airfoil.Neutralpoints()
Airfoil.plotairfoil()
print(Airfoil.CalculateEI())
print(Airfoil.curvatures(30000, 0))

EI = 5129815431.075207
epsilon = 3224999*23/EI
print(epsilon)