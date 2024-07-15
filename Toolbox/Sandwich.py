import numpy as np


class Sandwich:
    def __init__(self, laminate1, laminate2, core, Loads, Strains):
        self.laminate1 = laminate1
        self.laminate2 = laminate2
        self.core = core

        self.Loads = Loads
        self.Strains = Strains
        self.h = self.laminate1.h + self.laminate2.h + self.core.h

        # ABD matrix assigned upon initialisation:
        self.CalculateEquivalentABD()
        self.CalculateEquivalentProperties()

    def CalculateEquivalentABD(self):
        L1ABD = self.laminate1.CalculateCoreABD(self.core.h)
        L2ABD = self.laminate2.CalculateCoreABD(self.core.h)

        totalABD = L1ABD + L2ABD

        self.ABD_matrix = totalABD
        return totalABD

    def CalculateEquivalentProperties(self):
        # Here we calculate the engineering constants (or equivalent properties):
        self.Ex = (self.ABD_matrix[0, 0] * self.ABD_matrix[1, 1] - self.ABD_matrix[0, 1] ** 2) / (
                    self.h * self.ABD_matrix[1, 1])
        self.Ey = (self.ABD_matrix[0, 0] * self.ABD_matrix[1, 1] - self.ABD_matrix[0, 1] ** 2) / (
                    self.h * self.ABD_matrix[0, 0])

        self.vxy = self.ABD_matrix[0, 1] / self.ABD_matrix[1, 1]
        self.vyx = self.ABD_matrix[0, 1] / self.ABD_matrix[0, 0]

        self.Gxy = self.ABD_matrix[2, 2] / self.h
        return

    def FailureAnalysis(self):
        # Check wrinkling and FPF:
        FPFFI = self.LaminateFPF()
        print('first ply failure FI:', FPFFI)

        return

    def LaminateFPF(self):
        # TODO: implement correct load distribution between face sheets!
        laminate_Loads = [0.5 * x for x in self.loads]
        self.laminate1.Loads = laminate_Loads
        maxFI = self.laminate1.FailureAnalysis()[2]

        self.laminate2.Loads = laminate_Loads
        maxFI = self.laminate2.FailureAnalysis()[2]

        return maxFI

    def SandwichWrinkingShear(self):
        t_face = self.laminate1.h
        ABD_matrix45 = self.laminate1.rotated_ABD(np.deg2rad(45))

        # In this case we check for shear buckling:
        vxy = ABD_matrix45[0, 1] / ABD_matrix45[1, 1]
        vyx = ABD_matrix45[0, 1] / ABD_matrix45[0, 0]
        D11f = ABD_matrix45[3, 3]

        # And we calculate the effective E modulus of the face sheet:
        E_face = (12 * (1 - vxy * vyx) * D11f) / (t_face ** 3)

        G_45 = self.core.Gxbarz(np.deg2rad(45))
        Ez = self.core.coreproperties['Ez']
        t_core = self.core.thickness

        # Check which formula to use for wrinkling:
        symthick = self.SymThickWrinkling(Ez, t_face, E_face, G_45)

        asymthin = self.SymThinWrinkling(Ez, t_core, t_face, E_face, G_45)

        if symthick < asymthin:
            Nwrinkle = symthick
            print("SymThickWrinkling")
        else:
            Nwrinkle = asymthin
        return Nwrinkle

    def SymThickWrinkling(self, Ez, t_face, E_f, G_45):
        Ns_w = 0.43 * t_face * (E_f * Ez * G_45) ** (1 / 3)
        return Ns_w

    def SymThinWrinkling(self, Ez, t_core, t_face, E_f, G_45):
        Ns_w = (0.816 * np.sqrt(E_f * Ez * t_face ** 3 / t_core) + G_45 * t_core / 6)
        return Ns_w

    def AsymThickWrinkling(self, Ez, t_core, t_face, E_f):
        Ns_w = 0.33 * t_face * E_f * np.sqrt(Ez * t_face / (E_f * t_core))
        return Ns_w

    def AsymThinWrinkling(self, Ez, t_core, t_face, E_f, G_45):
        Ns_w = 0.33 * t_face ** (3 / 2) * np.sqrt(Ez * E_f / t_core)
        return Ns_w


