import numpy as np

np.set_printoptions(linewidth=300, precision = 3)

class Lamina:
    def __init__(self, t, theta, E1, E2, G12, v12, z0=None, z1=None, Sigma = None, Epsilon = None):
        self.t = t
        self.E1 = E1  # Modulus of elasticity in material direction 1
        self.E2 = E2  # Modulus of elasticity in material direction 2
        self.G12 = G12  # Shear modulus in plane 12
        self.v12 = v12  # Poisson's ratio in plane 12
        self.z0 = None
        self.z1 = None
        self.Epsilon = Epsilon
        self.Sigma = Sigma

        theta_rad = np.radians(theta)
        m = np.cos(theta_rad)  # Cosine of the angle of the material direction 1 with the x-axis
        n = np.sin(theta_rad)  # Sine of the angle of the material direction 1 with the x-axis
        # Now let's translate the equations to Python using numpy operations

        v21 = v12 * E2 / E1

        Q = 1 - v12 * v21
        Q11 = E1 / Q
        Q22 = E2 / Q
        Q12 = v12 * E2 / Q
        Q66 = G12

        Qxx = Q11 * m ** 4 + 2 * (Q12 + 2 * Q66) * (m ** 2) * (n ** 2) + Q22 * (n ** 4)
        Qxy = (Q11 + Q22 - 4 * Q66) * (m ** 2) * (n ** 2) + Q12 * (m ** 4 + n ** 4)
        Qyy = Q11 * (n ** 4) + 2 * (Q12 + 2 * Q66) * (m ** 2) * (n ** 2) + Q22 * (m ** 4)
        Qxs = (Q11 - Q12 - 2 * Q66) * n * (m ** 3) + (Q12 - Q22 + 2 * Q66) * (n ** 3) * m
        Qys = (Q11 - Q12 - 2 * Q66) * m * (n ** 3) + (Q12 - Q22 + 2 * Q66) * (m ** 3) * n
        Qss = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * (n ** 2) * (m ** 2) + Q66 * (n ** 4 + m ** 4)

        self.Q = np.array([[Qxx, Qxy, Qxs],
                      [Qxy, Qyy, Qys],
                      [Qxs, Qys, Qss]])

        # This Q is the ROTATED (global) q matrix, we can invert it to get the compliance matrix

        self.S = np.linalg.inv(self.Q)

    # We use the following method to carry out stress analysis for the laminate:
    def StressAnalysis(self):
        self.Sigma = np.zeros((3, 1))
        self.Sigma = self.Q @ self.Epsilon
        return self.Sigma

    # Finally, we can carry out failure analysis to figure out whether a lamina has failed:
    def FailureAnalysis(self):
        pass




class Laminate:
    def __init__(self, laminas, Loads=None, Strains=None):
        self.laminas = laminas  # the layers have an order and a thickness, so find the thickness of laminate
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

        # Initialize A_ij as a zero matrix
        # Assuming Q is a 2D array, we need to know its size to initialize A_ij.
        # We're assuming all Q matrices are the same size, as an example we use 3x3.
        A_matrix = np.zeros((3, 3))
        B_matrix = np.zeros((3, 3))
        D_matrix = np.zeros((3, 3))

        #Per lamina we calculate the three matrices
        for lamina in self.laminas:
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
        if self.Loads is not None:
            self.Strains = np.zeros((6, 1))
            self.Strains = np.linalg.inv(self.ABD_matrix) @ self.Loads
        else:
            print('loads is nonetype')

    def CalculateLoads(self):
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

    # def CalculateLaminaStrains(self):
    #     self.CalculateStrains()
    #     Strains = self.Strains
    #     for i in self.laminas:
    #         e1 = Strains[0] - ((i.z0 + i.z1) / 2) * Strains[3]
    #         e2 = Strains[1] - ((i.z0 + i.z1) / 2) * Strains[4]
    #         e3 = Strains[2] - ((i.z0 + i.z1) / 2) * Strains[5]
    #         i.Epsilon = np.array([[e1], [e2], [e3]])

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

    # code the failure criteria per lamina

    # def PlotStrain(self)

    # def PlotStress(self)

#now we test the code:

s0 = Lamina(0.0002, 0, 150e9, 20e9, 5e9, 0.3)
s1 = Lamina(0.0002, 0, 150e9, 20e9, 5e9, 0.3)
s2 = Lamina(0.0002, 0, 150e9, 20e9, 5e9, 0.3)
s3 = Lamina(0.0002, 0, 150e9, 20e9, 5e9, 0.3)

# create the laminas list, for the laminate function:
laminas = [s0, s1, s2, s3]

# creating the laminate object:
laminate = Laminate(laminas)

# now we can apply loads to the laminate: (in Newton per meter or newtonmeter per meter)
laminate.Loads = np.array([[10000], [0], [0], [0], [0], [0]])

print(laminate.StressAnalysis())

