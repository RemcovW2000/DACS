import numpy as np
from tqdm import tqdm
class Laminate:
    def __init__(self, laminas, Loads=None, Strains=None):
        # laminas is a python list of lamina objects, which make up the laminate.
        # The laminate layers are ordened by the order of the list.
        self.laminas = laminas

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
        self.calculate_ABD()

    def calculate_ABD(self):
        # Initialize A_ij as a zero matrix

        # Initalizing the A, B and D matrix:
        A_matrix = np.zeros((3, 3))
        B_matrix = np.zeros((3, 3))
        D_matrix = np.zeros((3, 3))

        # Per lamina we calculate the three matrices
        for lamina in self.laminas:
            # First we recalculate the Q and S matrix of the lamina:
            lamina.calculate_QS()

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

    def calculate_strains(self):
        # First we RECALCULATE the ABD matrix -> this because based on the failurestate of the lamina,
        # They will have different Q matrices
        self.calculate_ABD()

        # Then we check if the loads are assigned and calculate the strains:
        if self.Loads is not None:
            self.Strains = np.zeros((6, 1))
            self.Strains = np.linalg.inv(self.ABD_matrix) @ self.Loads
        else:
            print('loads is nonetype')

        return self.Strains

    def get_strains(self):
        # This will give strains in the lamina at the 'current' state of the laminate for given loads
        # But will not change the attribute. This is a way to 'read' the strains with current failure state.
        Strains = np.linalg.inv(self.ABD_matrix) @ self.Loads
        return Strains

    def calculate_loads(self):
        # First we RECALCULATE the ABD matrix -> this because based on the failurestate of the lamina,
        # They will have different Q matrices
        self.calculate_ABD()

        # Then we check if the strains are assigned and calculate the strains:
        if self.Strains is not None:
            self.Loads = np.zeros((6, 1))
            self.Loads = self.ABD_matrix @ self.Strains
        else:
            print('loads is nonetype')

    #after calculating LAMINATE strains, we can find the strains per lamina:
    def calculate_lamina_strains(self):
        # To calcualte lamina strains we first need global strainsL
        self.calculate_strains()

        Strains = self.Strains

        # we go through all the lamina:
        for i in self.laminas:
            # Given the fact that strain is a linear gradient in the laminate, we can find
            # The max strain per lamina by finding picking the max between strain at z0 vs z1
            max1 = max(Strains[0] - i.z0 * Strains[3], Strains[0] - i.z1 * Strains[3], key=abs)
            max2 = max(Strains[1] - i.z0 * Strains[4], Strains[1] - i.z1 * Strains[4], key=abs)
            max3 = max(Strains[2] - i.z0 * Strains[5], Strains[2] - i.z1 * Strains[5], key=abs)
            i.Epsilon = np.array([max1, max2, max3])


    def calculate_equivalent_properties(self):
        # Here we calculate the engineering constants (or equivalent properties):
        Ex = (self.A_matrix[0, 0] * self.A_matrix[1, 1] - self.A_matrix[0, 1] ** 2) / (self.h * self.A_matrix[1, 1])
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
        return [Ex, Ey, vxy, vyx, Gxy], [E1b, E2b, G12b, v12b, v21b]

    def calculate_equivalent_properties_attheta(self, theta):

        # Here we calculate the engineering constants (or equivalent properties):
        Ex = (self.A_matrix[0, 0] * self.A_matrix[1, 1] - self.A_matrix[0, 1] ** 2) / (self.h * self.A_matrix[1, 1])
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
        return [Ex, Ey, vxy, vyx, Gxy], [E1b, E2b, G12b, v12b, v21b]

    def stress_analysis(self):
        # We need to make sure the lamina have strains:
        self.calculate_lamina_strains()

        # we need a method to store the stresses so we can check the stresses
        shape = (3, len(self.laminas))
        stresses = np.zeros(shape)
        for count, i in enumerate(self.laminas):
            # calling of the i.stress_analysis() method should also save the stresses as attributes
            stressesnonflat = i.stress_analysis()
            stressesflat = stressesnonflat.flatten()
            stresses[:, count] = stressesflat
        return stresses

    # carry out failure analysis for all lamina in laminate
    def failure_analysis(self):
        # We need to make sure the lamina have stresses:
        self.stress_analysis()

        # We make an array to track the failed lamina (which one failed):
        failedlamina = []

        # Initializing an array to save the failure factors:
        FailureFactors = []

        # We want to potentially save the lamina which failed, not useful in this assignment though.
        for count, lamina in enumerate(self.laminas):
            # Now run for the lamina, the failure analysis
            results = lamina.failure_analysis(lamina.Sigma)

            # If the failure of the lamina is 1 (for IFF) or 2 (for FF), the lamina has failed
            if results[0] >= 1:
                failedlamina.append(count)

            # set the correct index of the failurestate:
            self.FailureState[count] = lamina.FailureState
            FailureFactors.append(max(results[1], results[2]))

        # We save the maximum failure factor in any of the lamina, to calculate the next loadstep:
        maxfailurefactor = np.max(FailureFactors)

        return self.FailureState, failedlamina, maxfailurefactor

    def progressive_damage_analysis(self, loadingratio, loadincrement):
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
            FailureState, failedlamina, maxfailurefactor = self.failure_analysis()

            # If a lamina has failed, save these loads and strains
            if failedlamina:
                FailureLoadsList.append(Loads)
                FailureStrainsList.append(self.get_strains())

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

    def produce_failure_envelope(self, loadincrement):
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

            FailureLoads, FailureStrains = self.progressive_damage_analysis(loadingratio, loadincrement)

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
            self.reset_failure_state()
        return E22vsE12FPF, E22vsE12LPF, S22vsS12FPF, S22vsS12LPF, FailureStrainsList

    def reset_failure_state(self):
        # First we reset the failure state vector in the laminate:
        self.FailureState = np.zeros(len(self.laminas))

        # Then we also reset these in the lamina
        for lamina in self.laminas:
            lamina.FailureState = 0

        # Lastly, we recalculate the ABD matrix:
        self.calculate_ABD()
        return

