import numpy as np
from Laminate import Laminate
from Lamina import Lamina
import MP


class Stringer:
    def __init__(self, LaminateH = None, LaminateV = None, Width = None, Height = None, Position = np.array([0, 0])):\

        self.area = 0
        self.location = 0
        self.sigma = 0
        self.Ex = 0

        self.LaminateH = LaminateH
        self.LaminateV = LaminateV

        self.Width = Width
        self.Height = Height

        # # At initialisation we assign area and ex!
        # self.CalculateEAMembers()
        # self.CalculateA()
        # self.CalculateEx()

    def CalculateA(self):
        avertical = self.LaminateV.h * (self.Height - self.LaminateH.h)
        ahorizontal = self.LaminateH.h * self.Width
        self.area = avertical + ahorizontal
        return

    def CalculateEx(self):
        self.Ex = self.EAtotal/self.area
        return

    def CalculateEAMembers(self):
        # Vertical member EA:
        Eiv = 1/(self.LaminateV.h * self.LaminateV.ABD_matrix_inverse[0, 0])

        # print('Eiv:', Eiv)
        # The width here is the height of the stringer
        bv = self.Height
        tv = self.LaminateV.h
        EAv = Eiv * bv * tv

        # Horizontal member EA:
        Eih = 1/(self.LaminateH.h * self.LaminateH.ABD_matrix_inverse[0,0])
        # print('Eih:', Eih)
        # The width here is the width of the stringer
        bh = self.Width
        th = self.LaminateH.h
        EAh = Eih * bh * th

        self.EAv = EAv
        self.EAh = EAh
        self.EAtotal = EAv + EAh
        return EAv, EAh

    # We can use the EA equivalent properties of the vertical and horizontal members
    # In our axis system y = 0 is at the bottom of the t stringer (so where it connects to skin)
    def CalculateYBar(self):
        # Prerequisite is to find the EA of the members:
        self.CalculateEAMembers()
        # Use the formula from LECTURE 4 slide 7:
        teller = self.EAh * (0.5*self.LaminateH.h) + self.EAv * (0.5*self.Height + self.LaminateH.h)
        noemer = self.EAh + self.EAv
        ybar = teller / noemer
        self.ybar = ybar
        return ybar

    def FindEIEquivalent(self):
        # Prerequisite is to calculate neutral axis and the EAv and EAh
        self.CalculateYBar()

        # First do the vertical:
        Av = self.LaminateV.h * (self.Height-self.LaminateH.h)
        print('Area vertical member: ', Av)
        dv = abs((0.5 * self.Height + 0.5 * self.LaminateH.h) - self.ybar)

        Ev = self.EAv / Av
        EIv = Ev*((self.LaminateV.h * (self.Height-self.LaminateH.h)**3)/12 + Av * dv)
        print('EIv: ', EIv)

        # Then add the horizontal:
        # This formula is from lecture 4 slide
        E_bih = 12/(self.LaminateH.h * self.LaminateH.ABD_matrix_inverse[3, 3])
        di = abs(self.ybar - self.LaminateH.h * 0.5)
        EIh = E_bih * ((self.Width * self.LaminateH.h**3)/12 + self.LaminateH.h * self.Width * di**2)
        EIEquivalent = EIv + EIh

        self.EIEquivalent = EIEquivalent
        return EIEquivalent

    def FailureAnalysis(self):
        # 3 failure analyses:
        # 1. First ply failure:
        self.FPFFF = self.FPFanalysis()

        # 2. Global buckling:
        BucklingFF = self.BucklingAnalysis()

        # 3.
        # CripplingFFhorizontal =


    def FPFanalysis(self):
        # Depending on the EA of each member (horizontal or vertical) the load for each
        # Member is calculated:
        Nxh = self.Fx*(self.EAh/(self.EAv+self.EAh))/self.Width
        Nxv = self.Fx*(self.EAv/(self.EAv+self.EAh))/(self.Height - self.LaminateH.h)

        self.LaminateV.Loads = np.array([Nxv, 0, 0, 0, 0, 0])
        self.LaminateH.Loads = np.array([Nxh, 0, 0, 0, 0, 0])
        return


# V0 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
# V1 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
# V2 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
# V3 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
#
# H0 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
# H1 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
# H2 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
# H3 = Lamina(MP.t, 0, MP.elasticproperties, MP.failureproperties)
#
# # create the laminas list, for the laminate function:
# LaminasV = [V0, V1, V2, V3]
# LaminasH = [H0, H1, H2, H3]
#
# # creating the laminate object:
# LaminateH = Laminate(LaminasH)
# LaminateV = Laminate(LaminasV)
#
# TStringer1 = Stringer(LaminateH, LaminateV, 30, 30)