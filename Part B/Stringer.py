import numpy as np
from Laminate import Laminate
from Laminate import LaminateBuilder
from Lamina import Lamina
import MP


class Stringer:
    def __init__(self, LaminateH = None, LaminateV = None, Width = None, Height = None, Position = np.array([0, 0])):\

        self.area = 0
        self.location = 0
        self.sigma = 0
        self.Ex = 0
        self.Fx = 0
        self.FrameSpacing = 1000
        self.print = False
        self.rho = 0.00161 # placeholder! although this is the correct value
        self.Failure = False

        self.LaminateH = LaminateH
        self.LaminateV = LaminateV

        self.Width = Width
        self.Height = Height

        # At initialisation we assign area and ex!
        self.CalculateEAMembers()
        self.CalculateA()
        self.CalculateEx()
        self.FindEIEquivalent()
        self.Calculatempl()



    def CalculateA(self):
        avertical = self.LaminateV.h * (self.Height - self.LaminateH.h)
        ahorizontal = self.LaminateH.h * self.Width
        self.area = avertical + ahorizontal
        return

    def CalculateEx(self):
        self.Ex = self.EAtotal/self.area
        return

    def Calculatempl(self):
        self.mpl = self.rho*self.area
        if self.print == True:
            print(self.mpl)
        return self.mpl

    def CalculateEAMembers(self):
        # Vertical member EA:
        Ev = self.LaminateV.Ex
        # The width here is the height of the stringer
        bv = self.Height - self.LaminateV.h
        tv = self.LaminateV.h
        EAv = Ev * bv * tv

        # Horizontal member EA:
        Eh = self.LaminateH.Ex
        # The width here is the width of the stringer
        bh = self.Width
        th = self.LaminateH.h
        EAh = Eh * bh * th

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

        if self.print == True:
            print('height of the neutral axis: ', np.round(ybar, 2), 'mm')
        return ybar

    def FindEIEquivalent(self):
        # Prerequisite is to calculate neutral axis and the EAv and EAh
        self.CalculateYBar()

        # First do the vertical:
        Av = self.LaminateV.h * (self.Height-self.LaminateH.h)
        dv = abs((0.5 * self.Height + 0.5 * self.LaminateH.h) - self.ybar)

        Ev = self.EAv / Av
        EIv = Ev*((self.LaminateV.h * (self.Height-self.LaminateH.h)**3)/12 + Av * dv)

        # Then add the horizontal:
        # This formula is from lecture 4 slide
        E_bih = 12/(self.LaminateH.h * self.LaminateH.ABD_matrix_inverse[3, 3])
        di = abs(self.ybar - self.LaminateH.h * 0.5)
        EIh = E_bih * ((self.Width * self.LaminateH.h**3)/12 + self.LaminateH.h * self.Width * di**2)
        EIEquivalent = EIv + EIh

        self.EIEquivalent = EIEquivalent
        return EIEquivalent

    def FailureAnalysis(self):
        self.CalculateMemberLoads()
        # 3 failure analyses:
        # 1. First ply failure:
        self.FPFFIv, self.FPFFIh = self.FPFanalysis()
        # Returns 2 factors, one for horizontal, one for vertical laminate

        # 2. Global buckling:
        self.BucklingFI = self.BucklingAnalysis()

        # 3.
        self.CripplingFIv, self.CripplingFIh = self.CripplingAnalysis()

        FIs = [self.FPFFIv, self.FPFFIh, self.BucklingFI, self.CripplingFIv, self.CripplingFIh]
        if any(FI > 1 for FI in FIs):
            self.Failure = True
        return

    def BucklingAnalysis(self):

        # Column buckling of just the stringer: simply supported on either side
        L = self.FrameSpacing
        Pcr = np.pi**2 *self.EIEquivalent/(L**2)
        if self.Fx >= 0:
            FF = 0
        elif self.Fx < 0:
            FF = np.abs(self.Fx/Pcr)

        if self.print == True:
            print('buckling FF: ', np.round(FF, 2))
        return FF

    def MemberCripplingAnalysis(self, Laminate, Width):
        b = Width
        D_matrix = Laminate.D_matrix
        D11 = D_matrix[0, 0]
        D12 = D_matrix[0, 1]
        D22 = D_matrix[1, 1]
        D66 = D_matrix[2, 2]
        lamda = 5/12
        t = Laminate.h
        Nc = 12*D66/b**2

        Ncrit = Nc * ((b/t)**0.717)/1.63
        return Ncrit

    def CripplingAnalysis(self):
        # First horizontal:
        hwidth = (self.Width-self.LaminateV.h)/2
        Ncrith = self.MemberCripplingAnalysis(self.LaminateH, hwidth)
        if self.Nxh <= 0:
            FIhorizontal = np.abs(self.Nxh/Ncrith)
        elif self.Nxh > 0:
            FIhorizontal = 0

        vwidth = (self.Height-self.LaminateH.h)
        Ncritv = self.MemberCripplingAnalysis(self.LaminateV, vwidth)

        if self.Nxv <= 0:
            FIvertical = np.abs(self.Nxv/Ncritv)
        elif self.Nxv > 0:
            FIvertical = 0

        if self.print == True:
            print('FI for crippling, vertical, horizontal: ', np.round(FIvertical,2), np.round(FIhorizontal, 2))
        return FIvertical, FIhorizontal

    def CalculateMemberLoads(self):
        self.Nxh = self.Fx*(self.EAh/(self.EAv+self.EAh))/self.Width
        self.Nxv = self.Fx*(self.EAv/(self.EAv+self.EAh))/(self.Height - self.LaminateH.h)
        return
    def FPFanalysis(self):
        # Depending on the EA of each member (horizontal or vertical) the load for each
        # Member is calculated:
        self.Nxh = self.Fx*(self.EAh/(self.EAv+self.EAh))/self.Width
        self.Nxv = self.Fx*(self.EAv/(self.EAv+self.EAh))/(self.Height - self.LaminateH.h)

        self.LaminateV.Loads = np.array([self.Nxv, 0, 0, 0, 0, 0])
        self.LaminateH.Loads = np.array([self.Nxh, 0, 0, 0, 0, 0])

        FPFFFV = self.LaminateV.FailureAnalysis()[2]
        FPFFFH = self.LaminateH.FailureAnalysis()[2]

        if self.print == True:
            print('FI for FPF, vertical, horizontal: ', np.round(FPFFFV,2), np.round(FPFFFH, 2))
        return FPFFFV, FPFFFH


# ----------------------------------------------------------------------------------------------
# Stringer 1:
# ----------------------------------------------------------------------------------------------
# creating the laminate objects:
LaminateH1 = LaminateBuilder([45, -45, 0, 0, 0, 0, 0], True, True, 1)
LaminateV1 = LaminateBuilder([0, 0, 0, 0, 0, 45, -45], False, False, 1)

TStringer_1 = Stringer(LaminateH1, LaminateV1, 20, 20)


# ----------------------------------------------------------------------------------------------
# Stringer 2:
# ----------------------------------------------------------------------------------------------
# creating the laminate object:
LaminateH2 = LaminateBuilder([45, -45, 0, 0, 0, 0, 0], True, True, 1)
LaminateV2 = LaminateBuilder([0, 0, 0, 0, 0, 45, -45], False, False, 1)

TStringer_2 = Stringer(LaminateH2, LaminateV2, 20, 20)

# ----------------------------------------------------------------------------------------------
# Stringer 3:
# ----------------------------------------------------------------------------------------------
# creating the laminate object:
LaminateH3 = LaminateBuilder([45, -45, 0, 0, 0, 0, 0], True, True, 1)
LaminateV3 = LaminateBuilder([0, 0, 0, 0, 0, 45, -45], False, False, 1)

TStringer_3 = Stringer(LaminateH3, LaminateV3, 20, 20)