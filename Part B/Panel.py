import numpy as np
from Laminate import Laminate
from Lamina import Lamina

class Panel:
    def __init__(self, Laminate=None):
        self.start     = 0
        self.stop      = 0
        self.length    = 0
        self.B_1       = 0
        self.B_2       = 0
        self.sigma     = 0
        self.q_b       = 0
        self.Ns        = 0 
        self.Nx        = 0
        self.Ex        = 0 #NOTE: will be assigned according master skin
        self.Laminate  = Laminate 
        self.depth     = 0
        self.Failure   = False

    def FailureAnalysis(self):
        # first ply failure: 
        self.FPFFI = self.FPFanalysis()
        # global buckling:
        self.BucklingFI = self.BucklingAnalysis()

        # monolithic design
        self.MonolithicBucklingFI = self.BucklingAnalysis_MonoLithic()

        FIs = [self.FPFFI, self.BucklingFI, self.MonolithicBucklingFI]
        if any(FI > 1 for FI in FIs):
            self.Failure = True

    def FPFanalysis(self):
        # collecting loads
        self.Laminate.Loads = np.array([self.Nx, 0, self.Ns, 0, 0, 0])
        # FPF analysis
        FPF = self.Laminate.FailureAnalysis()[2]
        return FPF

    def BucklingAnalysis(self):
        # two solutions for Ncrit
        D11  = self.Laminate.D_matrix[0][0]
        D22  = self.Laminate.D_matrix[1][1]
        D12  = self.Laminate.D_matrix[0][1]
        D66  = self.Laminate.D_matrix[2][2]
        k = self.Ns/(self.Nx + 1e-20)
        a = self.depth
        b = self.length
        c = (a*k)/(b*np.pi**2)
        constant    = np.sqrt(9 + 65536/81 * c**2)
        numerator   = (np.pi**2)*(D11 + 2*(D12 + 2*D66)*(a/b)**2 + D22*(a/b)**4)
        denominator = (a**2)*(2 - 8192/81 * c**2)
        N0          = numerator/(denominator + 1e-20) * (5 - constant)
        BucklingFI  = abs(self.Nx/N0)
        return BucklingFI


    def BucklingAnalysis_MonoLithic(self):
        # two solutions for Ncrit
        D11  = self.Laminate.D_matrix[0][0]
        D22  = self.Laminate.D_matrix[1][1]
        D12  = self.Laminate.D_matrix[0][1]
        D66  = self.Laminate.D_matrix[2][2]
        k = self.Ns/(self.Nx + 1e-20)
        a = self.depth
        b = 1e30
        c = (a*k)/(b*np.pi**2)
        constant    = np.sqrt(9 + 65536/81 * c**2)
        numerator   = (np.pi**2)*(D11 + 2*(D12 + 2*D66)*(a/b)**2 + D22*(a/b)**4)
        denominator = (a**2)*(2 - 8192/81 * c**2)
        N0          = numerator/(denominator + 1e-20) * (5 - constant)
        BucklingFI  = abs(self.Nx/N0)
        return BucklingFI