import numpy as np
class Laminate:
    def __init__(self, laminas, Loads=None, Strains=None):
        self.laminas = laminas  # the layers have an order and a thickness, so find the thickness of laminate
        self.FailureState = np.zeros(len(self.laminas))
        # calculate total thickness
        # find layer start and end height
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
        self.CalculateABD()

    def CalculateABD(self):
        # Initialize A_ij as a zero matrix
        # Assuming Q is a 2D array, we need to know its size to initialize A_ij.
        # We're assuming all Q matrices are the same size, as an example we use 3x3.
        A_matrix = np.zeros((3, 3))
        B_matrix = np.zeros((3, 3))
        D_matrix = np.zeros((3, 3))

        #Per lamina we calculate the three matrices
        for lamina in self.laminas:
            # First we recalculate the Q and S matrix of the lamina:
            lamina.CalculateQS()

            # Calculate the difference (Z_k - Z_k-1)
            delta_Z = lamina.z1 - lamina.z0
            # Update A_ij by adding the product of Q(k) and the difference in Z
            A_matrix += lamina.Q * delta_Z

            delta_Z_squared = lamina.z1 ** 2 - lamina.z0 ** 2
            B_matrix += 1 / 2 * (lamina.Q * delta_Z_squared)

            delta_Z_cubed = lamina.z1 ** 3 - lamina.z0 ** 3
            D_matrix += 1 / 3 * (lamina.Q * delta_Z_cubed)

        self.A_matrix = A_matrix
        self.B_matrix = B_matrix
        self.D_matrix = D_matrix

        self.ABD_matrix = np.block([
            [A_matrix, B_matrix],
            [B_matrix, D_matrix]
        ])

    def CalculateStrains(self):
        # First we RECALCULATE the ABD matrix -> this because based on the failurestate of the lamina,
        # They will have different Q matrices
        self.CalculateABD()

        if self.Loads is not None:
            self.Strains = np.zeros((6, 1))
            self.Strains = np.linalg.inv(self.ABD_matrix) @ self.Loads
        else:
            print('loads is nonetype')

    def CalculateLoads(self):
        # First we RECALCULATE the ABD matrix -> this because based on the failurestate of the lamina,
        # They will have different Q matrices
        self.CalculateABD()

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
        for i in self.laminas:
            # Given the fact that strain is a linear gradient in the laminate, we can find
            # The max strain per lamina by finding picking the max between strain at z0 vs z1
            max1 = max(Strains[0] - i.z0 * Strains[3], Strains[0] - i.z1 * Strains[3], key=abs)
            max2 = max(Strains[1] - i.z0 * Strains[4], Strains[1] - i.z1 * Strains[4], key=abs)
            max3 = max(Strains[2] - i.z0 * Strains[5], Strains[2] - i.z1 * Strains[5], key=abs)
            i.Epsilon = np.array([max1, max2, max3])

    def CalculateEquivalentProperties(self):
        Ex = (self.A_matrix[0, 0] * self.A_matrix[1, 1] - self.A_matrix[0, 1] ** 2) / (self.h * self.A_matrix[1, 1])
        Ey = (self.A_matrix[0, 0] * self.A_matrix[1, 1] - self.A_matrix[0, 1] ** 2) / (self.h * self.A_matrix[0, 0])

        vxy = self.A_matrix[0, 1] / self.A_matrix[1, 1]
        vyx = self.A_matrix[0, 1] / self.A_matrix[0, 0]

        Gxy = self.A_matrix[2, 2] / self.h

        return [Ex, Ey, vxy, vyx, Gxy]

    # implement stress analysis per lamina -> calculate stresses per lamina with Q matrix
    # Since we have laminate epsilon
    def StressAnalysis(self):
        # We need to make sure the lamina have strains:
        self.CalculateLaminaStrains()

        # we need a method to store the stresses so we can check the stresses (maybe print them)
        shape = (3, len(self.laminas))
        stresses = np.zeros(shape)
        for count, i in enumerate(self.laminas):
            # calling of the i.Stressanalysis() method should also save the stresses as attributes
            stressesnonflat = i.StressAnalysis()
            stressesflat = stressesnonflat.flatten()
            stresses[:, count] = stressesflat
        return stresses

    #carry out failure analysis for all lamina in laminate
    def FailureAnalysis(self):
        # We need to make sure the lamina have stresses:
        self.StressAnalysis()
        for count, lamina in enumerate(self.laminas):
            # Now run for the lamina, the failure analysis
            lamina.FailureAnalysis(lamina.Sigma)
            self.FailureState[count] = lamina.FailureState
        return self.FailureState

    def ProgressiveDamageAnalysis(self, loadingratio, loadincrement):
        # We want to expose the laminate to a loading ratio, given as a numpy array
        # for example: np.array([[800], [70], [0], [0], [0], [0]])
        # We then want to normalize the list and apply loads in the ratio of the list
        normalized_loadingratio = loadingratio / np.max(np.abs(loadingratio))

        # We can save the point of FPF and LPF for later plotting
        LPF = False
        FPF = False

        # When failure happens we make a datapoint containing:
        # The lamina stresses
        # The global loads
        # The global strains
        # We do this at the lamina level!

        # Now we loop thru the load increments untill we reach last ply failure
        n = 0
        while LPF == False:
            # First we calculate the load for this iteration:
            Loads = normalized_loadingratio * n * loadincrement
            self.Loads = Loads
            # Then we run the actual failure analysis for the laminate, which also returns the failure state of the laminate:
            FailureState = self.FailureAnalysis()

            # We also want the first ply failure load:
            # And to know which ply fails!
            if FPF == False and np.any(FailureState == 1):
                #First ply failure is true now
                FPF = True
                FPF_load = Loads

            # Then we check whether we achieved full failure of all the lamina:
            if np.all(FailureState >= 1):
                LPF = True
                LPF_load = Loads

            n += 1
        return FPF_load, LPF_load

