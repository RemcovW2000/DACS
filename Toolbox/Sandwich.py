import numpy as np
from Toolbox.Laminate import Laminate
from Toolbox.Laminate import LaminateBuilder
from Toolbox.Core import Core
import Data.Panels
import copy

class Sandwich:
    def __init__(self, laminate1, laminate2, core, Loads = None, Strains = None):
        self.laminate1 = laminate1 # bottom laminate
        self.laminate2 = laminate2 # top laminate
        self.core = core
        self.sandwich = True

        self.Loads = Loads
        self.Strains = Strains
        self.h = self.laminate1.h + self.laminate2.h + self.core.h

        # ABD matrix assigned upon initialisation:
        self.CalculateEquivalentABD()
        self.CalculateEquivalentProperties()

    def CalculateEquivalentABD(self):
        '''
        Calculates the ABD matrix for the sandwich structure
        :return:
        '''
        # TODO: ADD TRANSVERSE SHEAR EFFECTS!!
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
        '''
        this is the general faulure analysis function

        Output should be a single scalar value, containing the failure indicator
        FI>0, when FI = 1 failure occurs at the applied load
        Loads must be assigned BEFORE calling this function!
        '''
        # First calculate loads on each facesheet:
        # This function also ASSIGNS LOADS TO THE LAMINATES!
        self.FaceSheetLoadDistribution()

        # then we can check first ply failure for both facesheets:
        FPFFI = self.LaminateFPF()
        print('first ply failure FI:', FPFFI)

        # Now we can also check facesheet wrinkling for the specific loadcase
        # NOTE: both facesheets must be analyzed seperately
        # for both facesheets find the principal directions and stresses:
        L1stresses, L1directions = self.laminate1.principal_stresses_and_directions()

        L2stresses, L2directions = self.laminate2.principal_stresses_and_directions()

        # using the directions we can find the material strength in that direction
        wrinklingFI1 = self.WrinklingAnalysis(L1stresses, L1directions, self.laminate1)
        wrinklingFI2 = self.WrinklingAnalysis(L2stresses, L2directions, self.laminate2)
        # Now use known formulas to find the failure indicators:
        print('WrinklingFIs: ', wrinklingFI1, wrinklingFI2)
        # Return largest FI!
        FImax = max([wrinklingFI1, wrinklingFI2, FPFFI])
        return FImax

    def LaminateFPF(self):
        # loads are assigned
        maxFI1 = self.laminate1.FailureAnalysis()
        maxFI2 = self.laminate2.FailureAnalysis()
        return max(maxFI1, maxFI2)

    def BucklingSchalingFactor(self, Ncrit):
        # Assuming k = 1 for now:
        k = 1
        Nxcrit = Ncrit/(1+Ncrit/(self.core.h*self.core.G))
        return Nxcrit

    def WrinklingAnalysis(self, laminatestresses, laminatedirections, laminate):
        # We obtain stresses and directions:
        # TODO: take into account assymetric laminates
        negative_values = laminatestresses[laminatestresses < 0]
        if len(negative_values) > 0:
            # Find the index of the negative value with the largest absolute value
            max_negative_index = np.where(laminatestresses == negative_values[np.argmax(np.abs(negative_values))])[0][0]
            direction = laminatedirections[max_negative_index]

            # now calculate E at the given angle:
            vx, vy = direction[0], direction[1]
            theta = np.arctan2(vy, vx)  # Calculate the angle theta in radians

            cos_theta = np.cos(theta)
            sin_theta = np.sin(theta)

            # Calculate 1/E_theta using the formula
            E_theta_inv = (
                    (cos_theta ** 4) / laminate.Ex +
                    (sin_theta ** 4) / laminate.Ey +
                    (2 * sin_theta ** 2 * cos_theta ** 2) / laminate.Gxy +
                    (2 * laminate.vxy * cos_theta ** 2 * sin_theta ** 2) / laminate.Ex
            )

            # Inverse to get E_theta
            E_theta = 1 / E_theta_inv

            # use the E_theta to find the FI in the given direction
            G_core = self.core.Gxbarz(np.deg2rad(theta))
            Ez = self.core.coreproperties['Ez']
            t_core = self.core.h
            t_face = laminate.h


            symthick = abs(self.SymThickWrinkling(Ez, t_face, E_theta, G_core))

            asymthin = abs(self.SymThinWrinkling(Ez, t_core, t_face, E_theta, G_core))
            if symthick < asymthin:
                Nwrinkle = symthick
            else:
                Nwrinkle = asymthin

            FI = abs(laminatestresses[max_negative_index]/Nwrinkle)

        else:
            FI = 0
        return FI

    def FaceSheetLoadDistribution(self):
        # Normal loads are as follows:
        Nx = self.Loads[0]
        Ny = self.Loads[1]

        # Divide normal loads between facesheets based on EA of facesheets

        # Assign Nx:
        Ext1 = self.laminate1.Ex*self.laminate1.h
        Ext2 = self.laminate2.Ex*self.laminate2.h

        Nx1 = Nx*(Ext1/(Ext1+Ext2))
        Nx2 = Nx*(Ext2/(Ext1+Ext2))

        # Assign Ny:
        Eyt1 = self.laminate1.Ey*self.laminate1.h
        Eyt2 = self.laminate2.Ey*self.laminate2.h

        Ny1 = Ny*(Eyt1/(Eyt1 + Eyt2))
        Ny2 = Ny*(Eyt2/(Eyt1 + Eyt2))

        Ns = self.Loads[2]

        # Divide shear loads between facesheets based on shear stifness GA of facesheets
        Gt1 = self.laminate1.Gxy*self.laminate1.h
        Gt2 = self.laminate2.Gxy*self.laminate2.h

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

if __name__ == '__main__':
    sandwich = copy.deepcopy(Data.Panels.Sandwiches['PanelWingRoot'])
    sandwich.Loads = [-10, 0, 2, 0, 0, 0]
    sandwich.FailureAnalysis()

