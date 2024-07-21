import copy

from Data.Airfoils import airfoilcoords
from scipy.interpolate import interp1d
from Toolbox.Laminate import LaminateBuilder
import matplotlib.pyplot as plt
from Toolbox.DamagedRegion import *
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

        self.topmembers = topmembers
        self.botmembers = botmembers
        self.members = topmembers + botmembers

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
                Dl = np.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)
                EA = Dl * member.panel.h * member.panel.Ex
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
        EIxxt, EIyyt, EIxyt = self.EI(self.topmembers, 100, 'top')
        EIxxb, EIyyb, EIxyb = self.EI(self.botmembers, 100, 'bot')
        self.EIxx = EIxxt + EIxxb
        self.EIyy = EIyyt + EIyyb
        self.EIxy = EIxyt + EIxyb
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
                    p2 = [xlist[i + 1]- self.xbar, self.Top_height_at_neutralax(xlist[i + 1]) if side == 'top' else self.Bottom_height_at_neutralax(xlist[i + 1])]

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
        e = -kx*(y-self.ybar) + ky*(x-self.xbar)
        return e

    def GenerateBooms(self, memberlist, side, Mx, My, npoints):
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

    def ShearFlow(self, memberlist, Sx, Sy, Brxr, Bryr, side):
        for member in memberlist if side == 'top' else reversed(memberlist):
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

    def ShearMomentSegments(self, segmentlist):
        # shear moment must be zero around any point (we'll take neutral point)
        moment = 0
        for segment in segmentlist:
            p1 = np.array(segment.p1)
            p2 = np.array(segment.p2)
            center = p1 + (p2 - p1)/2

            rx, ry = center[0], center[1]

            # force vector:
            F = (segment.qs * (p2-p1)/np.linalg.norm(p2-p1)) * np.linalg.norm(p2-p1)
            Fx, Fy = F[0], F[1]
            M = rx * Fy -ry * Fx
            moment += M
        return moment

    def ShearCorrection(self):
        moment = 0
        toppoints = []
        for member in self.topmembers:
            moment += self.ShearMomentSegments(member.segments)
            toppointsnew = [member.segments[i].p1 for i in range(len(member.segments))]
            toppoints += toppointsnew

        botpoints = []
        for member in self.botmembers:
            moment += self.ShearMomentSegments(member.segments)
            botpointsnew = [member.segments[i].p1 for i in range(len(member.segments))]
            botpoints += botpointsnew

        A_airfoil = self.calculate_area(toppoints, botpoints)
        qs0 = moment/(2*A_airfoil)

        for member in self.topmembers:
            for segment in member.segments:
                segment.qs = segment.qs + qs0

        for member in self.botmembers:
            for segment in member.segments:
                segment.qs = segment.qs + qs0

        return qs0

    def CalculateShearFlow(self, Sx, Sy, Mx, My):
        # First generate booms for top and bot:
        self.GenerateBooms(self.topmembers, 'top', Mx, My, 100)
        self.GenerateBooms(self.botmembers, 'bot', Mx, My, 100)
        Brxr_top, Bryr_top = self.ShearFlow(self.topmembers, Sx, Sy, 0, 0, 'top')
        self.ShearFlow(self.botmembers, Sx, Sy, Brxr_top, Bryr_top, 'bot')
        # now we've calculated all shear flows
        return

    def PlotShearFlow(self):
        allsegments = []
        for member in self.topmembers:
            newsegments = [(segment.p1, segment.p2, abs(segment.qs)) for segment in member.segments]
            allsegments.extend(newsegments)  # Use extend instead of append to flatten the list

        for member in self.botmembers:
            newsegments = [(segment.p1, segment.p2, abs(segment.qs)) for segment in member.segments]
            allsegments.extend(newsegments)  # Use extend instead of append to flatten the list

        # Extract values to normalize
        values = [value for _, _, value in allsegments]
        max_value = max(values)

        # Normalize values
        normalized_values = [value / max_value for value in values]

        # Create a color map
        cmap = plt.get_cmap('viridis')

        # Plot the lines
        fig, ax = plt.subplots()
        for (p1, p2), normalized_value in zip(
                [(p1, p2) for p1, p2, value in allsegments],
                normalized_values):
            x1, y1 = p1
            x2, y2 = p2
            ax.plot([x1, x2], [y1, y2], color=cmap(normalized_value), linewidth=2)

        # Set equal scaling
        ax.set_aspect('equal', adjustable='box')

        # Add a color bar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max_value))
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label='Value')

        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_title('Lines with Value-Based Coloring')
        plt.show()
        return


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




if __name__ == '__main__':

    Laminate = LaminateBuilder([0], False, True, 1)
    Member = Member(Laminate)
    sparlocations = [80, 200]
    topmembers = [copy.deepcopy(Member) for _ in range(len(sparlocations)+1)]
    botmembers = [copy.deepcopy(Member) for _ in range(len(sparlocations)+1)]

    type = 'e395-il'
    Airfoil = Airfoil(type, 1, 300, sparlocations, topmembers, botmembers)
    Airfoil.Neutralpoints()
    Airfoil.plotairfoil()
    print(Airfoil.CalculateEI())
    print(Airfoil.curvatures(30000, 0))
    Mx = 30000
    My = 0
    Sx = 0
    Sy = 100
    Airfoil.CalculateShearFlow(Sx, Sy, Mx, My)
    Airfoil.ShearCorrection()
    Airfoil.PlotShearFlow()