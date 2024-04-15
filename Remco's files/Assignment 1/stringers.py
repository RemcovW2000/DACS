# We'll start with a cross sectional analysis for a T stringer
import numpy as np
from Laminate import Laminate
from Lamina import Lamina


class Stringer:
    def __init__(self, Ex, LaminateH = None, LaminateV = None, Width = None, Height = None, Position = np.array([0, 0])):

        self.LaminateH = LaminateH
        self.LaminateV = LaminateV

        self.Width = Width
        self.Height = Height



        self.CalculateEAMembers()   # This calculates equivalent EA for all the members

        self.area = 0
        self.location = 0
        self.sigma = 0
        self.Ex = 0

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
