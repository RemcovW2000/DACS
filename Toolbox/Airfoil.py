import copy

import numpy as np
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from Data.Airfoils import airfoilcoords
from scipy.interpolate import interp1d
from Toolbox.Laminate import LaminateBuilder
import matplotlib.pyplot as plt
from Toolbox.DamagedRegion import *
from Toolbox.Member import Member
from matplotlib.colors import Normalize, TwoSlopeNorm
from scipy.linalg import lstsq
import math
from Toolbox.section import section
from Toolbox.helperclasses import boom, segment

'''
Shear flow has to be fixed,

current hypothesis is that the shear moments are not being calculated correctly!

segments p1 and p2 are defined

the sign of qs may not by definition coincide with the direction of the segment (from p1 to p2)
'''

class Airfoil:
    def __init__(self, airfoilname, thickness, chordlength, sparlocations, topmembers, botmembers, sparmembers, trpanel = None, trstart = None, trend = None, brpanel = None, brstart = None, brend = None):
        """
        We make the assumption that the laminate does not change within a member!
        For now we'll assume the whole structure is sandwich panel?

        :param sparlocations: locations from LE to TE of spars
        :param airfoilname:
        :param thickness:
        :param chordlength:
        :param Cl:
        :param moments:
        - list of moments -> there are 2, around both the x and y axis
        :param shearforces:
        - list of shear forces on the profile -> there are again 2, in the direction of both axes

        :param trpanel: top reinforcement panel. This is the full panel of the reinforced area.
        :param brpanel: bottom reinforcement panel, again this is the full panel of the reinforced area
        :params trstart, trend, brstart, brend: top- and bottom reinforcement start and en points, in mm from the leading
        edge
        """
        self.print = False
        self.plot = False
        self.airfoilname = airfoilname
        self.thickness = thickness
        self.chordlength = chordlength

        self.topmembers = topmembers
        self.botmembers = botmembers
        self.members = topmembers + botmembers

        # make coordinates:
        self.topcoordinates_base = airfoilcoords[airfoilname]['top']
        self.bottomcoordinates_base = airfoilcoords[airfoilname]['bottom']
        self.top_bot_coordinates()

        self.sparlocations = sparlocations
        self.sparmembers = sparmembers

        # reinforcement data:
        self.trpanel = trpanel
        self.trstart = trstart
        self.trend = trend
        self.brpanel = brpanel
        self.brstart = brstart
        self.brend = brend

        self.segmentlength = 1 #length of segments in mm
        self.structuralidealisation = False
        self.Mx = 0
        self.My = 0
        self.Sx = 0
        self.Sy = 0
        self.yshear = 0
        self.xshear = 0

        # Initialize attributes that will be computed later
        self.xbar = None
        self.ybar = None
        self.EIxx = None
        self.EIyy = None
        self.EIxy = None

        # find the member locations upon initialisation: startcoord, endcoord, etc:
        self.FindMemberLocations()
        # assign the submembers to the members which have reinforcement:
        self.assignsubmembers()

# ------------------------------------------------------------------------------------------------------
# Helper functions, not called in the algorithms directly, but in functions

    def top_bot_coordinates(self):
        # transform the coordinates as neccesary:
        topcoordinates = self.topcoordinates_base
        botcoordinates = self.bottomcoordinates_base

        topcoordinates = [[element * self.chordlength for element in sublist] for sublist in topcoordinates]
        botcoordinates = [[element * self.chordlength for element in sublist] for sublist in botcoordinates]
        self.topcoordinates = topcoordinates
        self.botcoordinates = botcoordinates

        return self.topcoordinates, self.botcoordinates

    def calculate_area(self, top_points, bottom_points):
        """
        Calculate the area enclosed by the airfoil section, defined by the top and bottom surfaces.

        :param top_points: List of lists representing the points on the upper surface [(x1, y1), (x2, y2), ...]
        :param bottom_points: List of lists representing the points on the lower surface [(x1, y1), (x2, y2), ...]
        :return: The area enclosed by the box-shaped object.
        """

        # Ensure points are in correct order to form a closed polygon
        # Top points go from left to right and bottom points from right to left to form the perimeter
        polygon_points = top_points + bottom_points

        # Calculate area using the shoelace formula
        n = len(polygon_points)
        area = 0.0
        for i in range(n):
            x1, y1 = polygon_points[i]
            x2, y2 = polygon_points[(i + 1) % n]
            area += x1 * y2 - x2 * y1

        return abs(area) / 2.0

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

        for i, member in enumerate(self.sparmembers):
            member.startcoord = [self.sparlocations[i], self.Top_height_at(self.sparlocations[i])]
            member.endcoord = [self.sparlocations[i], self.Bot_height_at(self.sparlocations[i])]
        return

    def Top_height_at(self, x):
        xlist = [point[0] for point in self.topcoordinates]
        ylist = [point[1] for point in self.topcoordinates]
        top_interp = interp1d(xlist, ylist, kind='linear', bounds_error=False, fill_value="extrapolate")
        return float(top_interp(x))

    def Bot_height_at(self, x):
        xlist = [point[0] for point in self.botcoordinates]
        ylist = [point[1] for point in self.botcoordinates]
        bot_interp = interp1d(xlist, ylist, kind='linear', bounds_error=False, fill_value="extrapolate")
        return float(bot_interp(x))

    def Top_height_at_neutralax(self, x):
        return self.Top_height_at(x) - self.ybar

    def Bottom_height_at_neutralax(self, x):
        return self.Bot_height_at(x) - self.ybar

    def segment_I(self, t, P1, P2):
        """"
        Calculates Ixy, Ixx ,Iyy around neutral axis
        """
        # first Ixx:
        x1, y1 = P1
        x2, y2 = P2

        x = x1 + (x2-x1)/2
        y = y1 + (y2-y1)/2
        A = np.sqrt((x2-x1)**2 + (y2-y1)**2) * t
        if x1 == x2:
            h = abs(y2 - y1)
            Ixx = abs(0.125*t*h**3 + A * y **2)
            Iyy = abs(0.125*h*t**3 + A * x **2)
            Ixy = abs(A * x * y)
        elif y1 == y2:
            h = t
            b = abs(x2 - x1)
            Ixx = abs(0.125 * b * h ** 3 + A * y ** 2)
            Iyy = abs(0.125 * h * b ** 3 + A * x ** 2)
            Ixy = abs(A * x * y)
        else:
            b = (y1 - (x1 * y2 / x2)) / (1 - x1 / x2)
            a = (y2 - b) / x2
            Ixx = abs(t * np.sqrt(1 + a ** 2) * (
                        a ** 2 * (x2 ** 3 - x1 ** 3) / 3 + a * b * (x2 ** 2 - x1 ** 2) + b ** 2 * (x2 - x1)))

            # Ixy in a similar way:
            Ixy = t*np.sqrt(1+a**2) * ((x2**3 - x1**3)*a/3 + (x2**2-x1**2)*b/2)

            # We can simply switch x and y, and obtain Iyy using the same equations:
            y1, x1 = P1
            y2, x2 = P2

            if x2 == 0:
                x2 +=1e-15

            b = (y1 - (x1 * y2 / x2)) / (1 - x1 / x2)
            a = (y2 - b) / x2
            Iyy = abs(t * np.sqrt(1 + a ** 2) * (
                        a ** 2 * (x2 ** 3 - x1 ** 3) / 3 + a * b * (x2 ** 2 - x1 ** 2) + b ** 2 * (x2 - x1)))

        return Ixx, Iyy, Ixy

# ------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------
# Algorithm functions, called directly in the sequence for cross sectional analysis

    def assignsubmembers(self):
        """
        This function makes submembers and assigns them to the correct members:
        :return:
        """

        # -> check if there is even a reinforcement
        if self.trpanel:
            # print(self.trpanel)
            # print(self.trstart)
            # this means there is a top reinforcement!
            for member in self.topmembers:
                # for each member, see if the member contains a submember:
                # if the member contains reinforcement:
                if member.startcoord[0] < self.trstart < member.endcoord[0]:
                    # this means the start is in the member, is the end in?
                    member.submember = self.trpanel
                    if member.endcoord[0] > self.trend:
                        # this means the end is also in the member, assign member submember:
                        member.submember_start = self.trstart
                        member.submember_end = self.trend
                    else:
                        # end is not in member:
                        member.submember_start = self.trstart
                        member.submember_end = member.endcoord[0]

                # check if end is in member:
                elif member.startcoord[0] < self.trend < member.endcoord[0]:
                    # this means the end is in but the start is not
                    member.submember = self.trpanel
                    member.submember_start = member.startcoord[0]
                    member.submember_end = self.trend

        if self.brpanel:
            # this means there is a top reinforcement!
            for member in self.botmembers:
                # for each member, see if the member contains a submember:
                # if the member contains reinforcement:
                if member.startcoord[0] < self.brstart < member.endcoord[0]:
                    # this means the start is in the member, is the end in?
                    member.submember = self.brpanel
                    if member.endcoord[0] > self.brend:
                        # this means the end is also in the member, assign member submember:
                        member.submember_start = self.brstart
                        member.submember_end = self.brend
                    else:
                        # end is not in member:
                        member.submember_start = self.brstart
                        member.submember_end = member.endcoord[0]

                # check if end is in member:
                elif member.startcoord[0] < self.brend < member.endcoord[0]:
                    # this means the end is in but the start is not
                    member.submember = self.brpanel
                    member.submember_start = member.startcoord[0]
                    member.submember_end = self.brend
        if self.print:
            print('--------------------------------------------------------------------------------')
            print('Member reinforcements:')
            print('Topmembers:')
            for member in self.topmembers:
                print('start and end coordinates of member:')
                print(member.startcoord, member.endcoord)
                print('start and end x coordinates of reinforcement:')
                print(member.submember_start, member.submember_end)
                print('--------------------------------------------------------------------------------')
            print('Botmembers:')
            for member in self.botmembers:
                print('start and end coordinates of member:')
                print(member.startcoord, member.endcoord)
                print('start and end x coordinates of reinforcement:')
                print(member.submember_start, member.submember_end)
                print('--------------------------------------------------------------------------------')
        return

    def Neutralpoints(self):
        # find the neutral points -> x and y coordinate of neutral axes around x and y axes:
        EAytot = 0
        EAxtot = 0
        EAtot = 0
        for member in self.topmembers:
            npoints = 10  # TODO: placeholder!

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
                EA = Dl * member.h(x) * member.Ex(x)
                Dx = (p1[0] + p2[0]) / 2
                Dy = (p1[1] + p2[1]) / 2

                EAy = EA * Dy
                EAx = EA * Dx

                EAytot += EAy
                EAxtot += EAx
                EAtot += EA

        for member in self.botmembers:
            npoints = 10  # TODO: placeholder!

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
                Dl = np.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)
                EA = Dl * member.h(x) * member.Ex(x)
                Dx = (p1[0] + p2[0]) / 2
                Dy = (p1[1] + p2[1]) / 2

                EAy = EA * Dy
                EAx = EA * Dx

                EAytot += EAy
                EAxtot += EAx
                EAtot += EA

        for i, member in enumerate(self.sparmembers):
            # Spars are vertical:
            x = self.sparlocations[i]
            top = self.Top_height_at(x)
            bot = self.Bot_height_at(x)

            # avg height:
            y = (top + bot) / 2
            EA = member.panel.h * (abs(top - bot)) * member.panel.Ex # dont need to calculate, no reinforcement in spar!
            EAy = EA * y
            EAx = EA * x

            EAtot += EA
            EAxtot += EAx
            EAytot += EAy

        ybar = EAytot / EAtot
        xbar = EAxtot / EAtot

        self.xbar = xbar
        self.ybar = ybar
        # pass xbar and ybar to member:
        allmembers = self.topmembers + self.botmembers + self.sparmembers
        for member in allmembers:
            member.xbar = xbar
            member.ybar = ybar
        return xbar, ybar

    def plotairfoil(self):

        # Coloring inside of the airfoil:
        xlist = np.linspace(0.1, self.chordlength, 100)
        top = [self.Top_height_at(x) for x in xlist]
        bot = [self.Bot_height_at(x) for x in xlist]


        plt.figure(figsize=(10, 6))
        plt.fill_between(np.linspace(0, self.chordlength, len(xlist)), top, bot,
                         color='lightblue', alpha=0.5)
        def plotmembers(memberlist, color, label=None):
            for idx, member in enumerate(memberlist):
                xlist = [boom.location[0] + self.xbar for boom in member.booms]
                ylist = [boom.location[1] + self.ybar for boom in member.booms]

                # Place the text near the middle of each member
                mid_index = len(member.booms) // 2
                xlabel = xlist[mid_index] - 10
                ylabel = ylist[mid_index] + 5

                stacking_text = ', '.join(map(str, member.panel.stackingsequence))  # Converts list to a readable string
                plt.text(xlabel, ylabel, stacking_text, color=color, fontsize=5)
                plt.text(xlabel -5, ylabel, idx, color='black', fontsize=8)

                # Only set the label for the first member in the list
                plt.plot(xlist, ylist, label=label if idx == 0 else None, color=color)
            return

        # Plotting top members:
        plotmembers(self.topmembers, 'orange', label='Top Members')

        # Plotting bottom members:
        plotmembers(self.botmembers, 'blue', label='Bottom Members')

        # Plotting spar members:
        plotmembers(self.sparmembers, 'green', label='Spar Members')

        # plot reinforcement:
        if self.trpanel:
            trxlist = np.linspace(self.trstart, self.trend, 20)
            trylist = [self.Top_height_at(x) for x in trxlist]
            plt.plot(trxlist, trylist, label='Top reinforcement', color='black', linewidth=2)

        if self.brpanel:
            brxlist = np.linspace(self.brstart, self.brend, 20)
            brylist = [self.Bot_height_at(x) for x in brxlist]
            plt.plot(brxlist, brylist, label='Bottom reinforcement', color='black', linewidth=2)

        plt.axhline(y=self.ybar, linestyle='--', label = 'neutral bending axis')
        plt.axvline(x=self.xbar, linestyle='--', label = 'neutral bending axis')
        plt.axis('equal')
        plt.xlabel('X values')
        plt.ylabel('Heights')
        plt.title('Top and Bottom Heights vs. X values')
        plt.legend()
        plt.grid(True)
        return plt

    def CalculateEI(self):
        if self.xbar is None or self.ybar is None:
            raise ValueError(
                "Neutral points (xbar, ybar) must be calculated before calculating EI. Call Neutralpoints() first.")
        # for bot and top members we have a function:
        EIxxt, EIyyt, EIxyt = self.EI(self.topmembers, 100, 'top')
        EIxxb, EIyyb, EIxyb = self.EI(self.botmembers, 100, 'bot')

        # for the spars we do the following:
        EIxxs = 0
        EIyys = 0
        EIxys = 0
        for member in self.sparmembers:
            Ixx, Iyy, Ixy = self.segment_I(member.panel.h, member.startcoord, member.endcoord) # no reinforcement so no need for member.Ex(x)
            EIxx = Ixx * member.panel.Ex # no reinforcement in spar so no need for member.Ex(x)
            EIyy = Iyy * member.panel.Ex # no reinforcement in spar
            EIxy = Ixy * member.panel.Ex # no reinforcement in spar
            EIxxs += EIxx
            EIyys += EIyy
            EIxy += EIxy
        # Now add them all up to obtain the total:
        self.EIxx = EIxxt + EIxxb + EIxxs
        self.EIyy = EIyyt + EIyyb + EIyys
        self.EIxy = EIxyt + EIxyb + EIxys
        return EIxxt + EIxxb, EIyyt + EIyyb, EIxyt + EIxyb

    def EI(self, memberlist, npoints, side):
        EIxx, EIyy, EIxy = 0, 0, 0
        for member in memberlist:
            # in the topmembers list, we know the start and end coordinates:
            xlist = np.linspace(member.startcoord[0], member.endcoord[0], npoints)

            # now we need to find whether the member has a submember:
            if member.submember:
                # if so, we need the start and end points of the submember:

                start = member.submember_start + 1e-3  # add small value to have the point just inside the reinforcement

                # Find the index where x should be inserted
                indexstart = np.searchsorted(xlist, start)

                # Insert x into the array
                xlist = np.insert(xlist, indexstart, start)

                end = member.submember_end + 1e-3 # add small value to have the point just outside the reinforcement
                indexend = np.searchsorted(xlist, end)
                xlist = np.insert(xlist, indexend, end)
            else:
                pass

            for i, x in enumerate(xlist):
                if i == len(xlist) - 1:
                    break
                else:
                    # find the points with respect to neutral axes:
                    p1 = [x - self.xbar, self.Top_height_at_neutralax(x) if side == 'top' else self.Bottom_height_at_neutralax(x)]
                    p2 = [xlist[i + 1] - self.xbar, self.Top_height_at_neutralax(xlist[i + 1]) if side == 'top' else self.Bottom_height_at_neutralax(xlist[i + 1])]
                    # line segment starts at x and ends at x + 1

                Ilist = self.segment_I(member.h(x), p1, p2)
                EIxx += Ilist[0]*member.Ex(x)
                EIyy += Ilist[1]*member.Ex(x)
                EIxy += Ilist[2]*member.Ex(x)
        return EIxx, EIyy, EIxy

    def curvatures(self, Mx, My):
        EIxx = self.EIxx
        EIyy = self.EIyy
        EIxy = self.EIxy
        kx = (Mx * EIyy - My * EIxy) / (EIxx * EIyy - EIxy ** 2)
        ky = (My * EIxx - Mx * EIxy) / (EIxx * EIyy - EIxy ** 2)
        return kx, ky

    def AxialStrain(self, kx, ky, x, y):
        '''
        Calculates strain
        coordinates around the CENTROID/NEUTRAL POINT!!!
        '''
        e = -kx*(y) + ky*(x)
        return e

    def CalculateBoomAreas(self, p1, p2, member, kx, ky, x = None):
        # p1 and p2 are locations of booms with respect to neutral axes!
        # having 2 points, we obtain 2 strains and stresses
        e1 = self.AxialStrain(kx, ky, p1[0], p1[1])
        e2 = self.AxialStrain(kx, ky, p2[0], p2[1])

        # now find the stresses:
        if x is None:
            Ex = member.panel.Ex
        else:
            Ex = member.Ex(x)

        s1 = e1 * Ex
        s2 = e2 * Ex

        # based on this we can calculate areas for booms:
        b = np.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)
        if x is None:
            td = member.panel.h
        else:
            td = member.h(x)

        # IMPORTANT:
        # If the strains are both very close to 0 due to the points being close to the neutral axis, this calculation
        # will have numerical instabilities that deteriorate the accuracy of the solution. As a result, assume then
        # that the strain is equal in both booms (zero value thus equal), and assign the ratio between stresses as s1/s2 = 1:
        # check also 'Aircraft structures for engineering students' page 607.

        slow = 1e-11
        if e1 < slow or e2 < slow:
            B1 = (td * b / 6) * 3
            B2 = (td * b / 6) * 3
        else:
            B1 = (td * b / 6) * (2 + s2 / s1)
            B2 = (td * b / 6) * (2 + s1 / s2)

        # Calculate area ratio
        area_ratio = (B1 + B2) / (b * td)

        # Check bounds for the area ratio
        tolerance = 0.05  # Allowable deviation (5%)
        p1print = [p1[0]+self.xbar, p1[1] + self.ybar]
        p2print = [p2[0] + self.xbar, p2[1] + self.ybar]
        if not (1 - tolerance <= area_ratio <= 1 + tolerance):
            # TODO: fix this! when the magnitudes of strain are different due to higher or lower values for moment, the fixed value of 10e-11 does nothing!
            print(
                f"Warning: Area ratio out of bounds: {area_ratio:.3f} at point {p1print}, {p2print}, ybar = {self.ybar:.3f} "
                f"with member start {member.startcoord} and end {member.endcoord}. e1 = {e1}, e2 = {e2}"
            )
        return B1, B2, s1, s2

    def GenerateBooms(self, memberlist, side, Mx, My):
        '''
        Generate a bunch of skin segments with a thickness, start and end point, as well as sigma1, sigma2 and thickness
        '''
        kx, ky = self.curvatures(Mx, My)

        # looping through member lists and using previous function:
        for member in memberlist if side == 'top' else reversed(memberlist):
            start = np.array(member.startcoord)
            end = np.array(member.endcoord)

            if self.structuralidealisation:
                npoints = 2
            else:
                npoints = round(abs(np.linalg.norm(end-start)/self.segmentlength))

            # in the members lists, we know the start and end coordinates, make a list of x points along which we evaluate
            # the booms:
            xlist = np.linspace(member.startcoord[0] if side == 'top' else member.endcoord[0], member.endcoord[0] if side == 'top' else member.startcoord[0], npoints)


            if member.submember:
                # if so, we need the start and end points of the submember:

                if side == 'top':
                    start = member.submember_start + 1e-3  # add small value to have the point just inside the reinforcement
                    # Find the index where x should be inserted
                    indexstart = np.searchsorted(xlist, start)

                    # Insert x into the array
                    xlist = np.insert(xlist, indexstart, start)

                    end = member.submember_end + 1e-3 # add small value to have the point just outside the reinforcement
                    indexend = np.searchsorted(xlist, end)
                    xlist = np.insert(xlist, indexend, end)

                elif side == 'bot':
                    start = member.submember_start + 1e-3  # add small value to have the point just inside the reinforcement
                    # Find the index where x should be inserted for descending order
                    index = np.searchsorted(xlist[::-1], start)

                    # Insert x into the array
                    xlist = np.insert(xlist, len(xlist) - index, start)
                    end = member.submember_end + 1e-3  # add small value to have the point just outside the reinforcement

                    indexend = np.searchsorted(xlist[::-1], end)
                    xlist = np.insert(xlist, len(xlist) - indexend, end)
            else:
                pass

            boomplaceholder = boom()
            boomlist = [copy.deepcopy(boomplaceholder) for _ in xlist]

            for i, x in enumerate(xlist):
                # skip the last loop:
                if i == len(xlist) - 1:
                    break
                # then proceed:
                else:
                    # find p1 and p2:
                    p1 = [x-self.xbar, self.Top_height_at_neutralax(x) if side == 'top' else self.Bottom_height_at_neutralax(x)]

                    p2 = [xlist[i + 1]-self.xbar, self.Top_height_at_neutralax(xlist[i + 1]) if side == 'top' else self.Bottom_height_at_neutralax(xlist[i + 1])]


                    # feed appropriate x value to boom areas function:
                    B1, B2, s1, s2 = self.CalculateBoomAreas(p1, p2, member, kx, ky, x)
                    boomlist[i].B += B1
                    boomlist[i].Ex = member.Ex(x)
                    boomlist[i].location = p1
                    boomlist[i].Sigmax = s1
                    boomlist[i + 1].B += B2
                    boomlist[i + 1].Ex = member.Ex(x)
                    boomlist[i + 1].location = p2
                    boomlist[i + 1].Sigmax = s2
            member.booms = boomlist
            # TODO: write test that tests the number of booms and segments created
        return

    def GenerateSparBooms(self):
        '''
        Function generates booms for the spar members in the memberlists.
        :return:
        '''
        Mx = self.Mx
        My = self.My

        # Error check for moments, must be non zero
        if Mx == 0 and My == 0:
            raise ValueError("Moments Mx and My cannot both be zero. Analysis requires at least one non-zero moment.")

        kx, ky = self.curvatures(Mx, My)

        for i, member in enumerate(self.sparmembers):
            # apply the boom areas function to each segment in the member:
            start = np.array(member.startcoord)
            end = np.array(member.endcoord)
            if self.structuralidealisation:
                npoints = 2
            else:
                npoints = round(abs(np.linalg.norm(end - start) / self.segmentlength))
            ylocations = np.linspace(member.startcoord[1], member.endcoord[1], npoints)
            x = member.startcoord[0]

            boomplaceholder = boom()
            boomlist = [copy.deepcopy(boomplaceholder) for _ in ylocations]
            for i in range(len(ylocations)):
                # skip the last iteration in loop:
                if i == len(ylocations) - 1:
                    break
                # then proceed:
                # Find the coordinates of the start and end of a segment:
                p1 = [x - self.xbar, ylocations[i] - self.ybar]
                p2 = [x - self.xbar, ylocations[i + 1] - self.ybar]

                # Then find the boom areas and normal stresses:
                B1, B2, s1, s2 = self.CalculateBoomAreas(p1, p2, member, kx, ky)
                # Now add these to the boom objects in the list
                boomlist[i].B += B1
                boomlist[i].Ex = member.panel.Ex
                boomlist[i].location = p1
                boomlist[i].Sigmax = s1
                boomlist[i + 1].B += B2
                boomlist[i + 1].Ex = member.panel.Ex
                boomlist[i + 1].location = p2
                boomlist[i + 1].Sigmax = s2

            # Finally, assign the booms to the member:
            member.booms = boomlist
        return

    def Calculate_SectionShearFlow(self, memberlist, Sx, Sy):
        """ for one closed section, shear flow is calculated, and assigned to the members"""
        Brxr = 0
        Bryr = 0
        for member in memberlist:
            segmentexample = segment()
            segmentlist = [copy.deepcopy(segmentexample) for _ in range(len(member.booms) - 1)]

            # we need the 2nd moments of area first:

            EIxx = self.EIxx
            EIyy = self.EIyy
            EIxy = self.EIxy
            Term1 = -(Sx*EIxx - Sy*EIxy)/(EIxx*EIyy - EIxy**2)
            Term2 = -(Sy*EIyy - Sx*EIxy)/(EIxx*EIyy - EIxy**2)
            for i, boom in enumerate(member.booms):
                if i == len(member.booms)-1:
                    break
                else:
                    Brxr += boom.Ex * boom.B*boom.location[0]
                    Bryr += boom.Ex * boom.B*boom.location[1]
                    segmentlist[i].qs = Term1*Brxr + Term2*Bryr
                    segmentlist[i].p1, segmentlist[i].p2 = boom.location, member.booms[i+1].location
            member.segments = segmentlist
        return Brxr, Bryr

    def SectionShearFlows(self, Sx, Sy):
        self.GenerateBooms(self.topmembers, 'top', self.Mx, self.My)
        self.GenerateBooms(self.botmembers, 'bot', self.Mx, self.My)
        self.GenerateSparBooms()
        # initialize section list:
        self.sectionlist = []
        # There are as many sections as there are members in the top or bottom skin:
        for i, member in enumerate(self.topmembers):
            # exceptions for first and last sections:
            topmember = copy.deepcopy(self.topmembers[i])
            botmember = copy.deepcopy(self.botmembers[i])
            if i == 0:
                memberlist = [topmember, copy.deepcopy(self.sparmembers[i]), botmember]

            elif i == len(self.topmembers)-1:
                memberlist = [topmember, botmember,
                              copy.deepcopy(self.sparmembers[i-1])]

            # now the 'main loop'
            else:
                # NOTE: We make copies of all the members!
                memberlist = [topmember, copy.deepcopy(self.sparmembers[i]),
                              botmember, copy.deepcopy(self.sparmembers[i-1])]

            # now calculate the shear flow:
            # We then pass this member list and calculate the shear flow for the open cross section:
            self.Calculate_SectionShearFlow(memberlist, Sx, Sy)
            toppoints = [topmember.booms[i].location for i in range(len(topmember.booms))]
            botpoints = [botmember.booms[i].location for i in range(len(botmember.booms))]
            A = self.calculate_area(toppoints, botpoints)
            # determine if the current section is at a leading or trailing edge:
            leadingedge = False
            trailingedge = False
            if i == 0:
                leadingedge = True
            elif i == len(self.topmembers)-1:
                trailingedge = True

            # store in a section object for ease of access:
            currentsection = section(memberlist, Sx, Sy, leadingedge, trailingedge, A)
            self.sectionlist.append(currentsection)


        return

    def SectionShearCorrection(self):
        # condition matrix by using differnt units! currenly units are all N and mm
        # we can arbitrarily change the unit to match the magnitude of the numbers.

        # find the order of magnitude for the force intensities:
        sectionfortest = self.sectionlist[1]
        forceintesity = sectionfortest.prefactor * sectionfortest.deltawhole
        order_forceintensity = np.floor(np.log10(abs(forceintesity)))
        # print('order force intensity: ',order_forceintensity)

        # now find the order of magnitude for the areas:
        areasection = sectionfortest.A
        order_areasection = np.floor(np.log10(abs(areasection)))
        # print('order areasection: ',order_areasection)

        deltaorder_magnitude = abs(order_areasection - order_forceintensity)
        deltaorder_unit = np.round(deltaorder_magnitude/3)

        # figure out scaling factors:
        Ascalingfactor = 10**(-2*deltaorder_unit)
        qscalingfactor = 10**(deltaorder_unit)
        Mscalingfactor = 10**(-deltaorder_unit)
        radscalingfactor = 10**(order_forceintensity + deltaorder_unit)

        # now time to make the array:
        array = np.zeros((len(self.sectionlist) + 1, len(self.sectionlist) + 1))

        for i, section in enumerate(self.sectionlist):
            # print('section number:', i)
            if i == 0:
                print(section.deltaback)
                # special case: first section
                array[i, i] = section.prefactor * section.deltawhole * qscalingfactor
                array[i, i + 1] = -section.prefactor * section.deltaback * qscalingfactor
            elif i == len(self.sectionlist)-1:
                # print(section.deltafront)
                # special case: last section
                array[i, i] = section.prefactor * section.deltawhole * qscalingfactor
                array[i, i-1] = -section.prefactor * section.deltafront * qscalingfactor
            else:
                # print(section.deltafront)
                # print(section.deltaback)
                # all other sections:
                array[i, i] = section.prefactor * section.deltawhole * qscalingfactor
                array[i, i + 1] = -section.prefactor * section.deltaback * qscalingfactor
                array[i, i - 1] = -section.prefactor * section.deltafront * qscalingfactor

            array[i, len(self.sectionlist)] = -radscalingfactor

        # add last equation:
        array[-1, : -1] = [2*section.A * Ascalingfactor for section in self.sectionlist]
        # define vector:
        vector = np.zeros(len(self.sectionlist) + 1)
        for i, section in enumerate(self.sectionlist):
            vector[i] = -section.int_qb_ds_over_t * section.prefactor * qscalingfactor
        vector[-1] = sum([-section.ShearMoment * Mscalingfactor for section in self.sectionlist])
        
        # moment arms of shear forces wrt shear center (neutral point):
        y_arm = self.yshear - self.ybar
        x_arm = self.xshear - self.xbar
        self.Moment_applied = self.Sx*y_arm + -self.Sy * x_arm
        vector[-1] += self.Moment_applied * Mscalingfactor
        # print('--------------------------------------------------------------------------')
        # print('array:')
        # print(array)
        # print('vector:')
        # print(vector)
        # print('--------------------------------------------------------------------------')
        # we now solve the system
        x, residuals, rank, s = lstsq(array, vector)
        x[:-1] = x[:-1] / qscalingfactor
        x[-1] = x[-1] / radscalingfactor
        condnr = np.linalg.cond(array)
        # print('conditioning nr:', condnr)
        # print('solution vector:', x)
        # the shear flow corrections are as follows:
        correctionfactors = x[0:-1]
        # print('correction factors: ', correctionfactors)
        for i, Section in enumerate(self.sectionlist):
            Section.qs0 = correctionfactors[i]
            Section.ShearCorrection()

        if self.plot:
            for Section in self.sectionlist:
                Section.PlotShearFlow()
        return x

    def ShearSuperPosition(self):
        ''' this function superimposes the shear stresses found in the sections seperately, now in the original members'''

        for i, section in enumerate(self.sectionlist):
            if section.leadingedge:
                topindex = 0
                botindex = 2

                # now shear stress in the rib members: superimposed
                sparmemberindexes = [i]
                memberindexforspar = [1]

            elif section.trailingedge:
                topindex = 0
                botindex = 1

                sparmemberindexes = [i-1]
                memberindexforspar = [2]

            else:
                topindex = 0
                botindex = 2

                sparmemberindexes = [i, i - 1]
                memberindexforspar = [1, 3]
            self.topmembers[i] = section.members[topindex]
            self.botmembers[i] = section.members[botindex]

            # for the spars, the member class instance does not have segments yet, if it's the first time accessing the
            # spar member, we must copy it.

            # now shear stress in the spar members: superimposed
            for n, index in enumerate(sparmemberindexes):
                # if the spar member does not have segments, it means it has not been replaced by an
                # object from the analysis for the section, thus it can be replaced:
                if not self.sparmembers[index].segments:
                    self.sparmembers[index] = section.members[memberindexforspar[n]]
                else:
                    for k, segment in enumerate(self.sparmembers[index].segments):
                        # we must subtract the result as the qs in the next cell is always opposite!
                        segment.qs += -section.members[memberindexforspar[n]].segments[k].qs
        return

    def PlotShearFlow(self):
        allsegments = []
        for member in self.topmembers:
            newsegments = [(segment.p1, segment.p2, segment.qs) for segment in member.segments]
            allsegments.extend(newsegments)  # Use extend instead of append to flatten the list

        for member in self.botmembers:
            newsegments = [(segment.p1, segment.p2, segment.qs) for segment in member.segments]
            allsegments.extend(newsegments)  # Use extend instead of append to flatten the list

        for member in self.sparmembers:
            newsegments = [(segment.p1, segment.p2, segment.qs) for segment in member.segments]
            allsegments.extend(newsegments)

        # Extract values to normalize
        values = [value for _, _, value in allsegments]
        min_value = min(values)
        max_value = max(values)
        max_abs_value = max(abs(min_value), abs(max_value))

        # Determine the appropriate colormap and normalization
        if max_value <= 0:
            # All values are negative or zero, use absolute values for intensity
            cmap = plt.get_cmap('Blues')
            norm = Normalize(vmin=0, vmax=abs(min_value))
            normalized_values = [abs(value) for value in values]
        elif min_value >= 0:
            # All values are positive or zero
            cmap = plt.get_cmap('Reds')
            norm = Normalize(vmin=0, vmax=max_value)
            normalized_values = values
        else:
            # Values are both positive and negative, use symmetric normalization
            cmap = plt.get_cmap('coolwarm')
            norm = TwoSlopeNorm(vmin=-max_abs_value, vcenter=0, vmax=max_abs_value)
            normalized_values = values

        x_coords = []
        y_coords = []
        # Plot the lines
        fig, ax = plt.subplots()
        for (p1, p2), value in zip([(p1, p2) for p1, p2, value in allsegments], normalized_values):
            x1, y1 = p1
            x2, y2 = p2

            x_coords.extend([x1, x2])
            y_coords.extend([y1, y2])
            ax.plot([x1, x2], [y1, y2], color=cmap(norm(value)), linewidth=2)

        # Normalize the shear force vector to 0.25 * chord
        chord = self.chordlength
        vector_magnitude = (self.Sx ** 2 + self.Sy ** 2) ** 0.5
        normalized_sx = self.Sx / vector_magnitude * 0.1 * chord
        normalized_sy = self.Sy / vector_magnitude * 0.1 * chord

        # Add shear force vector using quiver
        ax.quiver(
            self.xshear - self.xbar, self.yshear - self.ybar, normalized_sx, normalized_sy,
            angles='xy', scale_units='xy', scale=1, color='black', label='Shear Force'
        )
        vector_tip_x = self.xshear - self.xbar + normalized_sx
        vector_tip_y = self.yshear - self.ybar + normalized_sy
        x_coords.append(vector_tip_x)
        y_coords.append(vector_tip_y)
        # Set equal scaling
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(min(x_coords) - 0.1 * chord, max(x_coords) + 0.1 * chord)
        ax.set_ylim(min(y_coords) - 0.05 * chord, max(y_coords) + 0.05 * chord)
        # Add a color bar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label='Value')
        ax.grid()
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_title('Shear flow, clockwise is positive, from the origin')
        return plt

    def PlotNormalStrain(self):
        '''
        we have to use the points given by points list for both top and bottom
        :return:
        '''
        topcoordinates = [[self.topcoordinates[i][0] - self.xbar, self.topcoordinates[i][1] - self.ybar] for i in range(len(self.topcoordinates))]
        botcoordinates = [[self.botcoordinates[i][0] - self.xbar, self.botcoordinates[i][1] - self.ybar] for i in range(len(self.botcoordinates))]
        coordinates = topcoordinates + botcoordinates

        x_coords = [point[0] for point in coordinates]
        y_coords = [point[1] for point in coordinates]
        kx, ky = self.curvatures(self.Mx, self.My)

        # having 2 points, we obtain 2 strains and stresses

        strains = [self.AxialStrain(kx, ky, point[0], point[1])*1e6 for point in coordinates]

        plt.figure(figsize=(8, 6))
        sc = plt.scatter(x_coords, y_coords, c=strains, cmap='viridis', s=20)

        # Add a color bar to indicate the value range
        plt.colorbar(sc, label='microstrain')

        # Normalize the axes
        plt.axis('equal')

        # Add labels and title
        plt.xlabel('X Coordinate from neutral point')
        plt.ylabel('Y Coordinate from neutral point')
        plt.title('Normal strain')

        # Show plot
        plt.show()
        return

    def PlotNormalForceIntensity(self):
        allsegments = []
        for member in self.topmembers:
            newsegments = [(segment.p1, segment.p2, member.booms[i].Sigmax*member.h(member.booms[i].location[0]-self.xbar)) for i, segment in enumerate(member.segments)]
            allsegments.extend(newsegments)  # Use extend instead of append to flatten the list

        for member in self.botmembers:
            newsegments = [(segment.p1, segment.p2, member.booms[i].Sigmax*member.h(member.booms[i].location[0]-self.xbar)) for i, segment in enumerate(member.segments)]
            allsegments.extend(newsegments)  # Use extend instead of append to flatten the list

        for member in self.sparmembers:
            newsegments = [(segment.p1, segment.p2, member.booms[i].Sigmax*member.panel.h) for i, segment in enumerate(member.segments)]
            allsegments.extend(newsegments)

        # Extract values to normalize
        values = [value for _, _, value in allsegments]
        min_value = min(values)
        max_value = max(values)
        max_abs_value = max(abs(min_value), abs(max_value))

        # Determine the appropriate colormap and normalization
        if max_value <= 0:
            # All values are negative or zero, use absolute values for intensity
            cmap = plt.get_cmap('Blues')
            norm = Normalize(vmin=0, vmax=abs(min_value))
            normalized_values = [abs(value) for value in values]
        elif min_value >= 0:
            # All values are positive or zero
            cmap = plt.get_cmap('Reds')
            norm = Normalize(vmin=0, vmax=max_value)
            normalized_values = values
        else:
            # Values are both positive and negative, use symmetric normalization
            cmap = plt.get_cmap('coolwarm')
            norm = TwoSlopeNorm(vmin=-max_abs_value, vcenter=0, vmax=max_abs_value)
            normalized_values = values

        # Plot the lines
        fig, ax = plt.subplots()
        for (p1, p2), value in zip([(p1, p2) for p1, p2, value in allsegments], normalized_values):
            x1, y1 = p1
            x2, y2 = p2
            ax.plot([x1, x2], [y1, y2], color=cmap(norm(value)), linewidth=2)

        # Set equal scaling
        ax.set_aspect('equal', adjustable='box')

        # Add a color bar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label='Nx')

        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_title('Normal force intensity (N/mm)')
        return plt
    def membercurvature(self, member, side):
        '''
        Assigns curvature to members. the reason we dont do this in the member is because the booms are located on linearly
        interpolated thus some form straight lines -> as a result, it would display curvature of 0.
        '''
        # find member start and end x coord:
        # Filter points based on the given x bounds
        lower_bound = member.startcoord[0]
        upper_bound = member.endcoord[0]

        # now filter the points in the correct list based on the lower and upper bound:
        if side == 'top':
            points = self.topcoordinates
        elif side == 'bot':
            points = self.botcoordinates
        else:
            raise ValueError('side not assigned, must be assigned either "top" or "bot"')
        filtered_points = [point for point in points if lower_bound <= point[0] <= upper_bound]

        # now find the radius:
        radius_max = 0

        for i in range (len(filtered_points)-2):
            p1 = filtered_points[i]
            p2 = filtered_points[i+1]
            p3 = filtered_points[i+2]
            radius = self.circle_radius_from_three_points(p1, p2, p3)
            if radius > radius_max:
                radius_max = radius
        member.R = radius_max
        return

    def AssignMemberCurvature(self):
        '''
        assigns curvature to all members in top and bottom:
        :return:
        '''
        for member in self.topmembers:
            self.membercurvature(member, side='top')

        for member in self.botmembers:
            self.membercurvature(member, side = 'bot')

        return

    def circle_radius_from_three_points(self, p1, p2, p3):
        # Extract coordinates
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3

        # Calculate the determinants
        D = np.linalg.det([[x1, y1, 1],
                           [x2, y2, 1],
                           [x3, y3, 1]])

        D_x = -np.linalg.det([[x1 ** 2 + y1 ** 2, y1, 1],
                              [x2 ** 2 + y2 ** 2, y2, 1],
                              [x3 ** 2 + y3 ** 2, y3, 1]])

        D_y = np.linalg.det([[x1 ** 2 + y1 ** 2, x1, 1],
                             [x2 ** 2 + y2 ** 2, x2, 1],
                             [x3 ** 2 + y3 ** 2, x3, 1]])

        C = -np.linalg.det([[x1 ** 2 + y1 ** 2, x1, y1],
                            [x2 ** 2 + y2 ** 2, x2, y2],
                            [x3 ** 2 + y3 ** 2, x3, y3]])

        # Calculate the center of the circle (a, b)
        a = -D_x / (2 * D)
        b = -D_y / (2 * D)

        # Calculate the radius R
        R = np.sqrt(a ** 2 + b ** 2 - C / D)

        return R
    def SolveStresses_CSA(self, moments, shearforces, center):
        '''
        Algorithm that solves for stresses, can be run after initialisation
        :return:
        '''
        self.Mx, self.My = moments
        self.Sx, self.Sy = shearforces
        # TODO: fix these dimensions so they are in the global (wing) FOR?
        self.xshear, self.yshear = center # relative to the local airfoil FOR

        # Functions for the normal stresses/deformations:
        self.Neutralpoints() # calculates neutral point
        self.CalculateEI() # then calculates EI around neutral point
        self.curvatures(self.Mx, self.My) # calculates curvatures (kx, ky)

        # Functions for shear stress solving
        self.SectionShearFlows(self.Sx, self.Sy)
        self.SectionShearCorrection()
        self.ShearSuperPosition()
        # we want to check the total shear force exerted by the section:
        self.ShearforceAnalysis()
        if self.plot:
            self.PlotShearFlow()
            self.PlotNormalForceIntensity()
        return

    def ShearforceAnalysis(self):
        Fxtop, Fytop, Mtop = self.force_moment_semgents(self.topmembers)
        Fxbot, Fybot, Mbot = self.force_moment_semgents(self.botmembers)
        Fxspar, Fyspar, Mspar = self.force_moment_semgents(self.sparmembers)

        Fx = Fxtop + Fxbot + Fxspar
        Fy = Fytop + Fybot + Fyspar
        M = Mtop + Mbot + Mspar
        return Fx, Fy, M

    def force_moment_semgents(self, memberlist):
        Fx = 0
        Fy = 0
        M = 0
        for member in memberlist:
            for segment in member.segments:
                # shear moments: (copied from section code)
                p1 = np.array(segment.p1)
                p2 = np.array(segment.p2)
                center = p1 + (p2 - p1) / 2
                rx, ry = center[0], center[1]

                # force vector:
                directionvector = p2 - p1
                # force by segments may not be consistent?!
                f = segment.qs * directionvector
                fx, fy = f[0], f[1]
                Fx += fx
                Fy += fy
                moment = ry * fx - rx * fy
                M += moment
        return Fx, Fy, M

    def FailureAnalysis_CSA(self):
        # first set the curvature of the members:
        self.AssignMemberCurvature()
        # then we can perform the failure analysis
        allmembers = self.topmembers + self.botmembers + self.sparmembers
        memberFIs = [member.CSA_FailureAnalysis() for member in allmembers]

        return max(memberFIs)



if __name__ == '__main__':

    Laminate = LaminateBuilder([45, -45, 0 ,90], True, True, 1)
    member = Member(Laminate)
    sparlocations = [80, 250]
    sparmembers = [copy.deepcopy(member) for _ in range(len(sparlocations))]
    topmembers = [copy.deepcopy(member) for _ in range(len(sparlocations)+1)]
    botmembers = [copy.deepcopy(member) for _ in range(len(sparlocations)+1)]

    reinforcementpaneltop = LaminateBuilder([45, -45, 0, 0, 0, 0, 0, 0, 90], True, True, 1)
    reinforcementpanelbot = LaminateBuilder([45, -45, 0, 0, 0, 0, 0, 0, 90], True, True, 1)

    topsubmember = Member(reinforcementpaneltop)
    botsubmember = Member(reinforcementpanelbot)
    reinforcement_start = 60
    reinforcement_end = 100

    type = 'NACA2410'
    Airfoil = Airfoil(type, 1, 300, sparlocations, topmembers, botmembers, sparmembers, topsubmember, 60, 100, botsubmember, 60, 100)
    Airfoil.SolveStresses_CSA([30000, 0], [0, 80], [300/4, 0])
    Airfoil.FailureAnalysis_CSA()

    # Assuming 'topshears' and 'xlist' are already defined from the loop
    topshears = []
    xlist = []
    for member in Airfoil.topmembers:
        for segment in member.segments:
            topshears.append(segment.qs)
            xlist.append((segment.p1[0] + segment.p2[0]) / 2)

    # Plot the shear forces
    plt.figure(figsize=(10, 6))
    plt.plot(xlist, topshears, label='Shear Force', linestyle='-')

    # Add labels, title, and legend
    plt.xlabel('x-coordinate (Chordwise Position)', fontsize=12)
    plt.ylabel('Shear Force (qs)', fontsize=12)
    plt.title('Chordwise Shear Force Distribution', fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)

    # Display the plot
    plt.show()
