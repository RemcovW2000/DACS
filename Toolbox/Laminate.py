import numpy as np
from tqdm import tqdm
from Toolbox.Lamina import Lamina
from Toolbox import MP
import copy
class Laminate:
    def __init__(self, laminas, Loads=None, Strains=None):
        # laminas is a python list of lamina objects, which make up the laminate.
        # The laminate layers are ordened by the order of the list, from bottom to top.
        self.laminas = laminas
        self.Loads   = Loads

        # The laminate also has a failure state, which is a list with the failure state of each lamina:
        self.FailureState = np.zeros(len(self.laminas))

        # We calculate the thickness and find layer start and end height:
        h = 0
        for i in laminas:
            # assign z0 and z1 for each layer:
            i.z0 = h
            i.z1 = h + i.t

            # keep track of total h
            h += i.t

        # now we subtract 0.5 times the height of the laminate to center it around z=0
        for i in laminas:
            i.z0 = i.z0 - 0.5 * h
            i.z1 = i.z1 - 0.5 * h
        self.h = h

        # We calculate the ABD matrix in initialisation
        self.CalculateABD()
        self.CalculateEquivalentProperties()

    def CalculateABD(self):
        # Initialize A_ij as a zero matrix

        # Initalizing the A, B and D matrix:
        A_matrix = np.zeros((3, 3))
        B_matrix = np.zeros((3, 3))
        D_matrix = np.zeros((3, 3))

        # Per lamina we calculate the three matrices
        for lamina in self.laminas:
            # First we recalculate the Q and S matrix of the lamina:
            lamina.CalculateQS()

            # Calculate the difference (Z_k - Z_k-1)
            delta_Z = lamina.z1 - lamina.z0
            # Update A_ij by adding the product of Q(k) and the difference in Z
            A_matrix += lamina.Q * delta_Z

            # Now the same for b and d matrices:
            delta_Z_squared = lamina.z1 ** 2 - lamina.z0 ** 2
            B_matrix += 1 / 2 * (lamina.Q * delta_Z_squared)

            delta_Z_cubed = lamina.z1 ** 3 - lamina.z0 ** 3
            D_matrix += 1 / 3 * (lamina.Q * delta_Z_cubed)

        # assign the matrices as attributes individually, this can be useful:
        # (but should be removed if this code should be used for high intensity applications)
        self.A_matrix = A_matrix
        self.B_matrix = B_matrix
        self.D_matrix = D_matrix

        # Save ABD matrix
        self.ABD_matrix = np.block([
            [A_matrix, B_matrix],
            [B_matrix, D_matrix]
        ])

        # We try the inversion but inversion is not so stable so we use an exception:
        try:
            self.ABD_matrix_inverse = np.linalg.inv(self.ABD_matrix)
        except:
            print('ABD may be singular')

    def CalculateStrains(self):
        # First we RECALCULATE the ABD matrix -> this because based on the failurestate of the lamina,
        # They will have different Q matrices
        self.CalculateABD()

        # Then we check if the loads are assigned and calculate the strains:
        if self.Loads is not None:
            self.Strains = np.zeros((6, 1))
            self.Strains = np.linalg.inv(self.ABD_matrix) @ self.Loads
        else:
            print('loads is nonetype')

        return self.Strains

    def GetStrains(self):
        # This will give strains in the lamina at the 'current' state of the laminate for given loads
        # But will not change the attribute. This is a way to 'read' the strains with current failure state.
        Strains = np.linalg.inv(self.ABD_matrix) @ self.Loads
        return Strains

    def CalculateLoads(self):
        # First we RECALCULATE the ABD matrix -> this because based on the failurestate of the lamina,
        # They will have different Q matrices
        self.CalculateABD()

        # Then we check if the strains are assigned and calculate the strains:
        if self.Strains is not None:
            self.Loads = np.zeros((6, 1))
            self.Loads = self.ABD_matrix @ self.Strains
        else:
            print('loads is nonetype')

    #after calculating LAMINATE strains, we can find the strains per lamina:
    def CalculateLaminaStrains(self):
        # To calcualte lamina strains we first need global strainsL
        self.CalculateStrains()

        Strains = self.Strains

        # we go through all the lamina:
        for i in self.laminas:
            # Given the fact that strain is a linear gradient in the laminate, we can find
            # The max strain per lamina by finding picking the max between strain at z0 vs z1
            max1 = max(Strains[0] - i.z0 * Strains[3], Strains[0] - i.z1 * Strains[3], key=abs)
            max2 = max(Strains[1] - i.z0 * Strains[4], Strains[1] - i.z1 * Strains[4], key=abs)
            max3 = max(Strains[2] - i.z0 * Strains[5], Strains[2] - i.z1 * Strains[5], key=abs)
            i.Epsilon = np.array([max1, max2, max3])


    def CalculateEquivalentProperties(self):
        # Here we calculate the engineering constants (or equivalent properties):
        self.Ex = (self.A_matrix[0, 0] * self.A_matrix[1, 1] - self.A_matrix[0, 1] ** 2) / (self.h * self.A_matrix[1, 1])
        Ey = (self.A_matrix[0, 0] * self.A_matrix[1, 1] - self.A_matrix[0, 1] ** 2) / (self.h * self.A_matrix[0, 0])

        vxy = self.A_matrix[0, 1] / self.A_matrix[1, 1]
        vyx = self.A_matrix[0, 1] / self.A_matrix[0, 0]

        Gxy = self.A_matrix[2, 2] / self.h

        D = np.linalg.inv(self.D_matrix)

        E1b = 12/(self.h**3 * D[0,0])
        E2b = 12/(self.h**3 * D[1,1])
        G12b = 12 / (self.h ** 3 * D[2,2])
        v12b = -D[0,1]/D[1,1]
        v21b = -D[0, 1] / D[0, 0]
        return [self.Ex, Ey, vxy, vyx, Gxy], [E1b, E2b, G12b, v12b, v21b]



    def StressAnalysis(self):
        # We need to make sure the lamina have strains:
        self.CalculateLaminaStrains()

        # we need a method to store the stresses so we can check the stresses
        shape = (3, len(self.laminas))
        stresses = np.zeros(shape)
        for count, i in enumerate(self.laminas):
            # calling of the i.Stressanalysis() method should also save the stresses as attributes
            stressesnonflat = i.StressAnalysis()
            stressesflat = stressesnonflat.flatten()
            stresses[:, count] = stressesflat
        return stresses

    # carry out failure analysis for all lamina in laminate
    def FailureAnalysis(self):
        # We need to make sure the lamina have stresses:
        self.StressAnalysis()

        # We make an array to track the failed lamina (which one failed):
        failedlamina = []

        # Initializing an array to save the failure factors:
        FailureFactors = []

        # We want to potentially save the lamina which failed, not useful in this assignment though.
        for count, lamina in enumerate(self.laminas):
            # Now run for the lamina, the failure analysis
            results = lamina.FailureAnalysis(lamina.Sigma)

            # If the failure of the lamina is 1 (for IFF) or 2 (for FF), the lamina has failed
            if results[0] >= 1:
                failedlamina.append(count)

            # set the correct index of the failurestate:
            self.FailureState[count] = lamina.FailureState
            FailureFactors.append(max(results[1], results[2]))

        # We save the maximum failure factor in any of the lamina, to calculate the next loadstep:
        maxfailurefactor = np.max(FailureFactors)
        return self.FailureState, failedlamina, maxfailurefactor

    def ProgressiveDamageAnalysis(self, loadingratio, loadincrement):
        # Normalize the loading ratio
        normalized_loadingratio = loadingratio / np.max(np.abs(loadingratio))

        # Last ply failure false at first:
        LPF = False

        # Initialize the failure loads and strains as empty lists
        FailureLoadsList = []
        FailureStrainsList = []

        n = 1
        while not LPF:
            # Calculate the load for this iteration
            Loads = normalized_loadingratio * n * loadincrement

            # Set the load attribute
            self.Loads = Loads

            # Run the failure analysis for the laminate with this new load
            FailureState, failedlamina, maxfailurefactor = self.FailureAnalysis()

            # If a lamina has failed, save these loads and strains
            if failedlamina:
                FailureLoadsList.append(Loads)
                FailureStrainsList.append(self.GetStrains())

            # Check whether full failure of all lamina has been achieved
            if np.all(FailureState >= 1):
                LPF = True

            if maxfailurefactor < 0.998:
                # The load should be increased based on the max failure factor observed:
                nnew = n*(1/maxfailurefactor)*0.999 + 1
                n = nnew
            else:
                n += 1

        # Convert lists to NumPy arrays for final output
        # If lists are empty, initialize arrays as (n,0) to avoid shape mismatch
        if FailureLoadsList:
            FailureLoads = np.hstack(FailureLoadsList)
            FailureStrains = np.hstack(FailureStrainsList)
        else:
            FailureLoads = np.empty((6, 0))
            FailureStrains = np.empty((6, 0))

        return FailureLoads, FailureStrains

    def ProduceFailureEnvelope(self, loadincrement):
        # We want to plot the stress and strain failure loads:
        angles = np.linspace(1, 360, 1440)
        E22vsE12FPF = []
        E22vsE12LPF = []

        S22vsS12FPF = []
        S22vsS12LPF = []

        FailureStrainsList = []

        for angle in tqdm(angles):
            loadingratio = np.array([[0],
                                     [np.cos(np.deg2rad(angle))],
                                     [np.sin(np.deg2rad(angle))],
                                     [0],
                                     [0],
                                     [0]])

            FailureLoads, FailureStrains = self.ProgressiveDamageAnalysis(loadingratio, loadincrement)

            # We save the individual points as tuples: This is for one load case:
            E22vsE12 = tuple(zip(FailureStrains[1], FailureStrains[2]))
            FailureStrainsList.append(E22vsE12)

            # now we save the FPF and LPF:
            E22vsE12FPF.append(E22vsE12[0])
            E22vsE12LPF.append(E22vsE12[-1])

            S22vsS12 = tuple(zip(FailureLoads[1], FailureLoads[2]))
            # print('S22 vs S12 failure loads:', S22vsS12, 'at angle:', angle)
            # Here we again take the FPF and LPF
            S22vsS12FPF.append(S22vsS12[0])
            S22vsS12LPF.append(S22vsS12[-1])
            self.ResetFailureState()
        return E22vsE12FPF, E22vsE12LPF, S22vsS12FPF, S22vsS12LPF, FailureStrainsList

    def ResetFailureState(self):
        # First we reset the failure state vector in the laminate:
        self.FailureState = np.zeros(len(self.laminas))

        # Then we also reset these in the lamina
        for lamina in self.laminas:
            lamina.FailureState = 0

        # Lastly, we recalculate the ABD matrix:
        self.CalculateABD()
        return

    def PrintAngles(self):
        angles = []
        for lamina in self.laminas:
            angles.append(lamina.theta)
        print(angles)

    def CalculateCoreABD(self, corethickness):
        # Initalizing the A, B and D matrix:
        A_matrix = np.zeros((3, 3))
        B_matrix = np.zeros((3, 3))
        D_matrix = np.zeros((3, 3))

        # Per lamina we calculate the three matrices
        for lamina in self.laminas:
            # First we recalculate the Q and S matrix of the lamina:
            lamina.CalculateQS()

            # Calculate the difference (Z_k - Z_k-1)
            z1 = lamina.z1 + corethickness/2 + self.h/2
            z0 = lamina.z0 + corethickness/2 + self.h/2
            delta_Z = z1 - z0
            # Update A_ij by adding the product of Q(k) and the difference in Z
            A_matrix += lamina.Q * delta_Z

            # Now the same for b and d matrices:
            delta_Z_squared = z1 ** 2 - z0 ** 2
            B_matrix += 1 / 2 * (lamina.Q * delta_Z_squared)

            delta_Z_cubed = z1 ** 3 - z0 ** 3
            D_matrix += 1 / 3 * (lamina.Q * delta_Z_cubed)

        # Save ABD matrix
        CoreABD = np.block([
            [A_matrix, B_matrix],
            [B_matrix, D_matrix]
        ])
        return CoreABD

    # Function to create the rotation matrix
    def rotation_matrix(self, theta):
        c = np.cos(theta)
        s = np.sin(theta)
        return np.array([[c ** 2, s ** 2, 2 * c * s],
                         [s ** 2, c ** 2, -2 * c * s],
                         [-c * s, c * s, c ** 2 - s ** 2]])

    # Function to transform the ABD matrix
    def rotated_ABD(self, theta):
        ABD = self.ABD_matrix
        T = self.rotation_matrix(theta)
        # Extending T to a 6x6 transformation matrix
        T_ext = np.zeros((6, 6))
        T_ext[:3, :3] = T
        T_ext[3:, 3:] = T

        # Transform the ABD matrix
        ABD_transformed = T_ext @ self.ABD_matrix @ T_ext.T
        return ABD_transformed

def LaminateBuilder(angleslist,symmetry, copycenter, multiplicity):
    if symmetry == True:
        if copycenter == True:
            angleslist = angleslist + angleslist[-1::-1]
        elif copycenter == False:
            angleslist = angleslist + angleslist[-2::-1]
    elif symmetry == False:
        angleslist = angleslist

    angleslist = angleslist * multiplicity

    # Define standard lamina:
    lamina = Lamina(MP.t, 45, MP.elasticproperties, MP.failureproperties)
    laminas = []

    # populate laminas list:
    for angle in angleslist:
        newlamina = copy.deepcopy(lamina)
        newlamina.theta = angle
        newlamina.CalculateQS()
        laminas.append(newlamina)

    laminate = Laminate(laminas)
    return laminate


