import numpy as np


class Sandwich:
    def __init__(self, laminate1, laminate2, core, Loads = None, Strains = None):
        self.laminate1 = laminate1 # bottom laminate
        self.laminate2 = laminate2 # top laminate
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
        # First calculate loads:
        self.FaceSheetLoadDistribution()

        # Check wrinkling and FPF:
        FPFFI = self.LaminateFPF()
        print('first ply failure FI:', FPFFI)

        return

    def LaminateFPF(self):
        # loads are assigned
        maxFI1 = self.laminate1.FailureAnalysis()[2]
        maxFI2 = self.laminate2.FailureAnalysis()[2]
        return max(maxFI1, maxFI2)

    def FaceSheetLoadDistribution(self):
        # Normal loads are as follows:
        Nx = self.Loads[0]
        Ny = self.Loads[1]

        # Divide normal loads between facesheets based on EA of facesheets

        # Assign Nx:
        Ext1 = self.laminate1.Ex*self.laminate1.t
        Ext2 = self.laminate2.Ex*self.laminate2.t

        Nx1 = Nx*(Ext1/(Ext1+Ext2))
        Nx2 = Nx*(Ext2/(Ext1+Ext2))

        # Assign Ny:
        Eyt1 = self.laminate1.Ey*self.laminate1.t
        Eyt2 = self.laminate2.Ey*self.laminate2.t

        Ny1 = Ny*(Eyt1/(Eyt1 + Eyt2))
        Ny2 = Ny*(Eyt2/(Eyt1 + Eyt2))

        Ns = self.Loads[2]

        # Divide shear loads between facesheets based on shear stifness GA of facesheets
        Gt1 = self.laminate1.Gxy*self.laminate1.t
        Gt2 = self.laminate2.Gxy*self.laminate2.t

        Ns1 = Ns*(Gt1/(Gt1 + Gt2))
        Ns2 = Ns*(Gt2/(Gt1 + Gt2))

        Mx = self.Loads[3]
        My = self.Loads[4]
        Ms = self.Loads[5]

        # moments are a bit more complicated but not strictly neccesary now
        # TODO: add facesheet loads due to moments

        self.laminate1.Loads = [Nx1, Ny1, Ns1, 0, 0, 0]
        self.laminate2.Loads = [Nx2, Ny2, Ns2, 0, 0, 0]
        return

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


