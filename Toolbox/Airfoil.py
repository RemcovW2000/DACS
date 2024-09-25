import copy

import numpy as np

from Data.Airfoils import airfoilcoords
from scipy.interpolate import interp1d
from Toolbox.Laminate import LaminateBuilder
import matplotlib.pyplot as plt
from Toolbox.DamagedRegion import *
from Toolbox.Member import Member
from scipy.integrate import simps
from matplotlib.colors import Normalize, TwoSlopeNorm
from scipy.linalg import lstsq
import math

class Airfoil:
    def __init__(self, airfoilname, thickness, chordlength, sparlocations, topmembers, botmembers, sparmembers):
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
        """
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

        self.FindMemberLocations()

        self.segmentlength = 1 #length of segments in mm
        self.structuralidealisation = True
        self.Mx = 0
        self.My = 0
        self.Sx = 0
        self.Sy = 0
        self.N0 = 0
        self.E0 = 0

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
            member.startcoord = [self.sparlocations[i], self.Bot_height_at(self.sparlocations[i])]
            member.endcoord = [self.sparlocations[i], self.Top_height_at(self.sparlocations[i])]
        return

    def Top_height_at(self, x):
        xlist = [point[0] for point in self.topcoordinates]
        ylist = [point[1] for point in self.topcoordinates]
        top_interp = interp1d(xlist, ylist, kind='linear', bounds_error=False, fill_value="extrapolate")
        return top_interp(x)

    def Bot_height_at(self, x):
        xlist = [point[0] for point in self.botcoordinates]
        ylist = [point[1] for point in self.botcoordinates]
        bot_interp = interp1d(xlist, ylist, kind='linear', bounds_error=False, fill_value="extrapolate")
        return bot_interp(x)

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
                EA = Dl * member.panel.h * member.panel.Ex
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
                EA = Dl * member.panel.h * member.panel.Ex
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
            EA = member.h * (abs(top - bot)) * member.panel.Ex
            EAy = EA * y
            EAx = EA * x

            EAtot += EA
            EAxtot += EAx
            EAytot += EAy

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
        plt.fill_between(np.linspace(0, self.chordlength, len(xlist)), top, bot,
                         color='lightblue', alpha=0.5)
        for spar in self.sparlocations:
            bot_at_rib = self.Bot_height_at(spar)
            top_at_rib = self.Top_height_at(spar)
            plt.plot([spar, spar], [bot_at_rib, top_at_rib], color='green', linewidth = 3, label='spar' if spar == self.sparlocations[0] else "")

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
        # for bot and top members we have a function:
        EIxxt, EIyyt, EIxyt = self.EI(self.topmembers, 100, 'top')
        EIxxb, EIyyb, EIxyb = self.EI(self.botmembers, 100, 'bot')
        # for the spars we do the following:
        EIxxs = 0
        EIyys = 0
        EIxys = 0
        for member in self.sparmembers:
            Ixx, Iyy, Ixy = self.segment_I(member.h, member.startcoord, member.endcoord)
            EIxx = Ixx * member.panel.Ex
            EIyy = Iyy * member.panel.Ex
            EIxy = Ixy * member.panel.Ex
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
            for i, x in enumerate(xlist):
                if i == len(xlist) - 1:
                    break
                else:
                    # find the points with neutral axes:
                    p1 = [x - self.xbar, self.Top_height_at_neutralax(x) if side == 'top' else self.Bottom_height_at_neutralax(x)]
                    p2 = [xlist[i + 1] - self.xbar, self.Top_height_at_neutralax(xlist[i + 1]) if side == 'top' else self.Bottom_height_at_neutralax(xlist[i + 1])]

                Ilist = self.segment_I(member.panel.h, p1, p2)
                EIxx += Ilist[0]*member.panel.Ex
                EIyy += Ilist[1]*member.panel.Ex
                EIxy += Ilist[2]*member.panel.Ex
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

    def GenerateBooms(self, memberlist, side, Mx, My):
        '''
        Generate a bunch of skin segments with a thickness, start and end point, as well as sigma1, sigma2 and thickness
        '''
        kx, ky = self.curvatures(Mx, My)
        def CalculateBoomAreas(p1, p2, member):
            # having 2 points, we obtain 2 strains and stresses
            e1 = self.AxialStrain(kx, ky, p1[0], p1[1])
            e2 = self.AxialStrain(kx, ky, p2[0], p2[1])

            # now find the stresses:
            Ex = member.panel.Ex
            s1 = e1 * Ex
            s2 = e2 * Ex

            # based on this we can calculate areas for booms:
            b = np.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)
            td = member.panel.h
            B1 = (td * b / 6) * (2 + s2 / s1)
            B2 = (td * b / 6) * (2 + s1 / s2)
            return B1, B2, s1, s2

        for member in memberlist if side == 'top' else reversed(memberlist):
            start = np.array(member.startcoord)
            end = np.array(member.endcoord)
            if self.structuralidealisation:
                npoints = 2
            else:
                npoints = round(abs(np.linalg.norm(end-start)/self.segmentlength))
            # in the topmembers list, we know the start and end coordinates:
            xlist = np.linspace(member.startcoord[0] if side == 'top' else member.endcoord[0], member.endcoord[0] if side == 'top' else member.startcoord[0], npoints)
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

                    B1, B2, s1, s2 = CalculateBoomAreas(p1, p2, member)
                    boomlist[i].B += B1
                    boomlist[i].location = p1
                    boomlist[i].Sigmax = s1
                    boomlist[i + 1].B += B2
                    boomlist[i + 1].location = p2
                    boomlist[i + 1].Sigmax = s2
            member.booms = boomlist
        return

    def GenerateSparBooms(self):
        Mx = self.Mx
        My = self.My
        kx, ky = self.curvatures(Mx, My)
        def CalculateBoomAreas(p1, p2, member):
            # having 2 points, we obtain 2 strains and stresses
            e1 = self.AxialStrain(kx, ky, p1[0], p1[1])
            e2 = self.AxialStrain(kx, ky, p2[0], p2[1])

            # now find the stresses:
            Ex = member.panel.Ex
            s1 = e1 * Ex
            s2 = e2 * Ex

            # based on this we can calculate areas for booms:
            b = np.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)
            td = member.panel.h
            B1 = (td * b / 6) * (2 + s2 / s1)
            B2 = (td * b / 6) * (2 + s1 / s2)
            return B1, B2, s1, s2

        for i, member in enumerate(self.sparmembers):
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
                # skip the last loop:
                if i == len(ylocations) - 1:
                    break
                # then proceed:
                # Find the coordinates of the stard and end of a segment:
                p1 = [x - self.xbar, ylocations[i] - self.ybar]
                p2 = [x - self.xbar, ylocations[i + 1] - self.ybar]

                # Then find the boom areas and normal stresses:
                B1, B2, s1, s2 = CalculateBoomAreas(p1, p2, member)

                # Now add these to the boom objects in the list
                boomlist[i].B += B1
                boomlist[i].location = p1
                boomlist[i].Sigmax = s1
                boomlist[i + 1].B += B2
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

            Ixx = self.EIxx/member.panel.Ex
            Iyy = self.EIyy/member.panel.Ex
            Ixy = self.EIxy/member.panel.Ex
            Term1 = -(Sx*Ixx - Sy*Ixy)/(Ixx*Iyy - Ixy**2)
            Term2 = -(Sy*Iyy - Sx*Ixy)/(Ixx*Iyy - Ixy**2)
            for i, boom in enumerate(member.booms):
                if i == len(member.booms)-1:
                    break
                else:
                    Brxr += boom.B*boom.location[0]
                    Bryr += boom.B*boom.location[1]
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
        for Section in self.sectionlist:
            Section.PlotShearFlow()
        # set up system of equations!

        # condition matrix by using differnt units! currenly units are all N and mm
        # we can arbitrarily change the unit to match the magnitude of the numbers.

        # find the order of magnitude for the force intensities:
        sectionfortest = self.sectionlist[1]
        forceintesity = sectionfortest.prefactor * sectionfortest.deltawhole
        order_forceintensity = np.floor(np.log10(abs(forceintesity)))
        print('order force intensity: ',order_forceintensity)

        # now find the order of magnitude for the areas:
        areasection = sectionfortest.A
        order_areasection = np.floor(np.log10(abs(areasection)))
        print('order areasection: ',order_areasection)

        deltaorder_magnitude = abs(order_areasection - order_forceintensity)
        deltaorder_unit = np.round(deltaorder_magnitude/3)

        # figure out scaling factors:
        Ascalingfactor = 10**(-2*deltaorder_unit)
        qscalingfactor = 10**(deltaorder_unit)
        Mscalingfactor = 10**(-deltaorder_unit)
        radscalingfactor = 10**(order_forceintensity + deltaorder_unit)

        array = np.zeros((len(self.sectionlist) + 1, len(self.sectionlist) + 1))

        for i, section in enumerate(self.sectionlist):
            if i == 0:
                array[i, i] = section.prefactor * section.deltawhole * qscalingfactor
                array[i, i + 1] = -section.prefactor * section.deltaback * qscalingfactor
            elif i == len(self.sectionlist)-1:
                array[i, i] = section.prefactor * section.deltawhole * qscalingfactor
                array[i, i-1] = -section.prefactor * section.deltafront * qscalingfactor
            else:
                array[i, i] = section.prefactor * section.deltawhole * qscalingfactor
                array[i, i + 1] = -section.prefactor * section.deltaback * qscalingfactor
                array[i, i - 1] = -section.prefactor * section.deltafront * qscalingfactor

            array[i, len(self.sectionlist)] = -1 * radscalingfactor

        # add last equation:
        array[-1, : -1] = [2*section.A * Ascalingfactor for section in self.sectionlist]
        # define vector:
        vector = np.zeros(len(self.sectionlist) + 1)
        for i, section in enumerate(self.sectionlist):
            vector[i] = -section.int_qb_ds_over_t * section.prefactor * qscalingfactor
        vector[-1] = sum([section.ShearMoment * Mscalingfactor for section in self.sectionlist])

        vector[-1] += (self.Sx*self.N0 - self.Sy * self.E0) * Mscalingfactor
        print('--------------------------------------------------------------------------')
        print('array:')
        print(array)
        print('vector:')
        print(vector)
        print('--------------------------------------------------------------------------')
        # we now solve the system
        x, residuals, rank, s = lstsq(array, vector)
        x[:-1] = x[:-1] / qscalingfactor
        x[-1] = x[-1] / radscalingfactor
        condnr = np.linalg.cond(array)
        print('conditioning nr:', condnr)
        print('solution vector:', x)
        # the shear flow corrections are as follows:
        correctionfactors = x[0:-1]
        print('correction factors: ', correctionfactors)
        for i, Section in enumerate(self.sectionlist):
            Section.qs0 = correctionfactors[i]
            Section.ShearCorrection()

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
        fig.colorbar(sm, ax=ax, label='Value')

        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_title('Shear flow, clockwise is positive, from the origin')
        plt.show()
        return

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
            newsegments = [(segment.p1, segment.p2, member.booms[i].Sigmax*member.panel.h) for i, segment in enumerate(member.segments)]
            allsegments.extend(newsegments)  # Use extend instead of append to flatten the list

        for member in self.botmembers:
            newsegments = [(segment.p1, segment.p2, member.booms[i].Sigmax*member.panel.h) for i, segment in enumerate(member.segments)]
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
        plt.show()
        return
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
            print('side not assigned, must be assigned either "top" or "bot"')
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
        performs sequence of
        :return:
        '''
        self.Mx, self.My = moments
        self.Sx, self.Sy = shearforces
        # TODO: fix these dimensions so they are in the global (wing) FOR?
        self.E0, self.N0 = center # now relative to neutral axis point/centroid

        # Functions for the normal stresses/deformations:
        self.Neutralpoints()
        self.CalculateEI()
        self.curvatures(self.Mx, self.My)

        # Functions for shear stress solving
        self.SectionShearFlows(self.Sx, self.Sy)
        self.SectionShearCorrection()
        self.ShearSuperPosition()
        # we want to check the total shear force exerted by the section:
        self.ShearforceAnalysis()
        self.PlotShearFlow()
        return

    def ShearforceAnalysis(self):
        Fx = 0
        Fy = 0
        for member in self.topmembers:
            for segment in member.segments:
                x1, y1 = segment.p1
                x2, y2 = segment.p2
                dx = x2 - x1
                dy = y2 - y1
                Fx += segment.qs * dx
                Fy += segment.qs * dy

        for member in self.botmembers:
            for segment in member.segments:
                x1, y1 = segment.p1
                x2, y2 = segment.p2
                dx = x2 - x1
                dy = y2 - y1
                print(dy)
                print('qs:', segment.qs)
                Fx += -segment.qs * dx
                Fy += -segment.qs * dy

        for member in self.sparmembers:
            for segment in member.segments:
                x1, y1 = segment.p1
                x2, y2 = segment.p2
                dx = x2 - x1
                dy = y2 - y1
                Fx += -segment.qs * dx
                Fy += -segment.qs * dy

        print('Fx from segments: ', Fx)
        print('Fy from segments: ', Fy)
        return

    def FailureAnalysis_CSA(self):
        # first set the curvature of the members:
        self.AssignMemberCurvature()
        # then we can perform the failure analysis
        allmembers = self.topmembers + self.botmembers + self.sparmembers
        memberFIs = [member.CSA_FailureAnalysis() for member in allmembers]

        return max(memberFIs)

class boom:
    def __init__(self):
        self.B = 0
        self.Ex = None
        self.Sigmax = None
        self.location = [None, None]
        self.qs = 0


class segment:
    def __init__(self):
        self.qs = None
        self.p1 = []
        self.p2 = []
class section:
    def __init__(self, members, Sx, Sy, leadingedge, trailingedge, A):
        self.members = members
        self.Sx, self.Sy = Sx, Sy
        self.leadingedge = leadingedge
        self.trailingedge = trailingedge
        self.A = A

        self.qs0 = None
        self.Corrected = False

        # Functions to run on initialisation:
        self.CalculateShearMoment()
        self.CalculateDeltas()
        self.Calculate_int_qb_ds_over_t()
        self.CalculateG()
        self.CalculatePrefactor()

    def CalculateShearMoment(self):
        moment = 0
        for member in self.members:
            moment += self.ShearMomentSegments(member.segments)
        self.ShearMoment = moment
        return self.ShearMoment

    def ShearMomentSegments(self, segmentlist):
        # shear moment must be zero around any point (we'll take neutral point)
        moment = 0
        for segment in segmentlist:
            p1 = np.array(segment.p1)
            p2 = np.array(segment.p2)
            center = p1 + (p2 - p1)/2
            rx, ry = center[0], center[1]

            # force vector:
            directionvector = p2-p1
            # force by segments may not be consistent?!

            F = segment.qs * directionvector
            Fx, Fy = F[0], F[1]
            M = rx * Fy -ry * Fx
            moment += M

        return moment

    def CalculateDeltas(self):
        if not self.leadingedge:
            frontmember = self.members[-1]

            difference = np.array(frontmember.endcoord) - np.array(frontmember.startcoord)

            arclength = np.linalg.norm(difference)
            self.deltafront = arclength / frontmember.panel.h
        if not self.trailingedge:
            backmember = self.members[1]
            difference = np.array(backmember.endcoord) - np.array(backmember.startcoord)

            arclength = np.linalg.norm(difference)
            self.deltaback = arclength / backmember.panel.h

        deltawhole = 0
        for member in self.members:
            # find arc length of member:
            segments = member.segments
            memberarclength = 0
            for segment in segments:
                # find the length of each segment and add to member arc length
                p1 = np.array(segment.p1)
                p2 = np.array(segment.p2)
                segmentlength = abs(np.linalg.norm(p2-p1))
                memberarclength += segmentlength

            member.arclength = memberarclength

            # add the delta of this member!
            deltawhole += memberarclength/member.panel.h

        self.deltawhole = deltawhole
        return

    def Calculate_int_qb_ds_over_t(self):
        productlist = []
        for member in self.members:
            productlistcurrent = [(segment.qs * np.linalg.norm(np.array(segment.p2)-np.array(segment.p1))) /member.panel.h for segment in member.segments]
            productlist += productlistcurrent

        integral = math.fsum(productlist)
        self.int_qb_ds_over_t = integral
        return self.int_qb_ds_over_t

    def CalculateG(self):
        # find the average g modulus of the cross section somehow?
        # TODO: find correct way to calculate G modulus
        # average wrt the arc length? -> weird, it should have something to do with it's contribution to the torsional stifness.
        Glist = [member.panel.Gxy * member.arclength for member in self.members]
        Arclengthlist = [member.arclength for member in self.members]
        self.G = sum(Glist)/sum(Arclengthlist)
        return self.G

    def CalculatePrefactor(self):
        self.prefactor = 1/(2 * self.G * self.A)
        return 1/(2 * self.G * self.A)

    def ShearCorrection(self):
        if not self.Corrected:
            for member in self.members:
                for segment in member.segments:
                    segment.qs += self.qs0
            self.Corrected = True
        else:
            print('WARNING, cross section has already been corrected! Correction not carried out.')
        return

    def PlotShearFlow(self):
        allsegments = []
        for member in self.members:
            newsegments = [(segment.p1, segment.p2, segment.qs) for segment in member.segments]
            allsegments.extend(newsegments)  # Use extend instead of append to flatten the list

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
        for (p1, p2), value in zip(
                [(p1, p2) for p1, p2, value in allsegments],
                normalized_values):
            x1, y1 = p1
            x2, y2 = p2
            ax.plot([x1, x2], [y1, y2], color=cmap(norm(value)), linewidth=2)

        # Set equal scaling
        ax.set_aspect('equal', adjustable='box')

        # Add a color bar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label='Value')

        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_title('Lines with Value-Based Coloring')
        plt.show()
        return




if __name__ == '__main__':

    Laminate = LaminateBuilder([45, -45, 0 ,90], True, True, 1)
    Member = Member(Laminate)
    sparlocations = [80, 200]
    sparmembers = [copy.deepcopy(Member) for _ in range(len(sparlocations))]
    topmembers = [copy.deepcopy(Member) for _ in range(len(sparlocations)+1)]
    botmembers = [copy.deepcopy(Member) for _ in range(len(sparlocations)+1)]

    type = 'NACA2410'
    Airfoil = Airfoil(type, 1, 300, sparlocations, topmembers, botmembers, sparmembers)
    Airfoil.SolveStresses_CSA([30000, 0], [0, 100], [-70, 0])
    # print(Airfoil.FailureAnalysis_CSA())
