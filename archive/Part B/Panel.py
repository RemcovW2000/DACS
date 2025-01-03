import numpy as np


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

    def failure_analysis(self):
        # first ply failure: 
        self.FPFFI = self.FPFanalysis()
        # global buckling:
        self.BucklingFI = self.buckling_analysis()

        # monolithic design
        self.MonolithicBucklingFI = self.buckling_analysis_MonoLithic()

        FIs = [self.FPFFI, self.BucklingFI, self.MonolithicBucklingFI]
        if any(FI > 1 for FI in FIs):
            self.Failure = True

    def FPFanalysis(self):
        # collecting loads
        self.Laminate.Loads = np.array([self.Nx, 0, self.Ns, 0, 0, 0])
        # FPF analysis
        FPF = self.Laminate.failure_analysis()[2]
        return FPF

    def buckling_analysis(self):
        # two solutions for Ncrit
        if self.Nx <= 1:
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
        else:
            BucklingFI = 0
        return BucklingFI


    def buckling_analysis_MonoLithic(self):
        # two solutions for Ncrit
        if self.Nx <= 1: 
            D11  = self.Laminate.D_matrix[0][0]
            D22  = self.Laminate.D_matrix[1][1]
            D12  = self.Laminate.D_matrix[0][1]
            D66  = self.Laminate.D_matrix[2][2]
            k = self.Ns/(self.Nx + 1e-20)
            a = self.depth
            b = 2.612e3
            c = (a*k)/(b*np.pi**2)
            constant    = np.sqrt(9 + 65536/81 * c**2)
            numerator   = (np.pi**2)*(D11 + 2*(D12 + 2*D66)*(a/b)**2 + D22*(a/b)**4)
            denominator = (a**2)*(2 - 8192/81 * c**2)
            N0          = numerator/(denominator + 1e-20) * (5 - constant)
            MonolithicBucklingFI  = abs(self.Nx/N0)
            #print('mono', MonolithicBucklingFI)
        else:
            MonolithicBucklingFI = 0
        # elif -0.1 < self.Nx < 0.1:
        #     D11  = self.Laminate.D_matrix[0][0]
        #     D22  = self.Laminate.D_matrix[1][1]
        #     D12  = self.Laminate.D_matrix[0][1]
        #     D66  = self.Laminate.D_matrix[2][2]
        #     k = self.Ns/(self.Nx + 1e-20)
        #     a = self.depth
        #     b = 0.1
        #     c = (a*k)/(b*np.pi**2)
        #     constant    = np.sqrt(9 + 65536/81 * c**2)
        #     numerator   = (np.pi**2)*(D11 + 2*(D12 + 2*D66)*(a/b)**2 + D22*(a/b)**4)
        #     denominator = (a**2)*(2 - 8192/81 * c**2)
        #     N0          = numerator/(denominator + 1e-20) * (5 - constant)
        #     MonolithicBucklingFI  = abs(self.Nx/N0)
        #     
        
        # else:
        #     MonolithicBucklingFI = 0

        return MonolithicBucklingFI 