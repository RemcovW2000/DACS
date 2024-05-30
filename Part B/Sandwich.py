import numpy as np
from tqdm import tqdm
from Lamina import Lamina
import MP
import copy
from Laminate import Laminate

class Core:
    def __init__(self, thickness, coreproperties):
        self.thickness = thickness
        self.coreproperties = coreproperties

    def Gxbarz(self, theta):
        Gxz = self.coreproperties['Gxz']
        Gyz = self.coreproperties['Gyz']
        sin_theta_squared = np.sin(theta) ** 2
        cos_theta_squared = np.cos(theta) ** 2
        Gxbarz = sin_theta_squared * Gyz + cos_theta_squared * Gxz
        return Gxbarz

    def Gybarz(self, theta):
        Gxz = self.coreproperties['Gxz']
        Gyz = self.coreproperties['Gyz']
        sin_theta_squared = np.sin(theta) ** 2
        cos_theta_squared = np.cos(theta) ** 2
        Gybarz = cos_theta_squared * Gyz + sin_theta_squared * Gxz
        return Gybarz

class Sandwich:
    def __init__(self, laminate1, laminate2, core, Loads):
        self.laminate1 = laminate1
        self.laminate2 = laminate2
        self.core = core

        self.Loads = [Loads[0], Loads[1], Loads[2], Loads[3], Loads[4], Loads[5]]

        # ABD matrix assigned upon initialisation:
        self.CalculateEquivalentABD()

    def CalculateEquivalentABD(self):
        L1ABD = self.laminate1.CalculateCoreABD(self.core.thickness)
        L2ABD = self.laminate2.CalculateCoreABD(self.core.thickness)

        totalABD = L1ABD + L2ABD

        self.ABD_matrix = totalABD
        return totalABD

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

class Member:
    def __init__(self, panel, loads, a, b):
        self.panel = panel # could be either laminate or sandwich, both should work
        self.loads = loads # total loads, not intensity!!!

        self.a = a # panel width
        self.b = b # panel depth

    def ShearBucklingFF(self):
        D11 = self.panel.ABD_matrix[3, 3]
        D12 = self.panel.ABD_matrix[3, 4]
        D22 = self.panel.ABD_matrix[4, 4]
        D66 = self.panel.ABD_matrix[5, 5]

        term1 = (9 * np.pi ** 4 * self.b) / (32 * self.a ** 3)
        term2 = (D11 + 2 * (D12 + 2 * D66)) * (self.a ** 2 / self.b ** 2)
        term3 = D22 * (self.a ** 4 / self.b ** 4)
        Nxypanel = 0.78 * term1 * (term2 + term3)
        return Nxypanel

    def PanelFF(self):
        # Calculate the failure state/factor of the panel itself, meaning FPF or crippling:
        # Assuming symmetric laminates, so the force is split between laminates:

        return None




