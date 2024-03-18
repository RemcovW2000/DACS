import numpy as np

np.set_printoptions(linewidth=300, precision = 3)

class Lamina:
    def __init__(self, t, theta, elasticproperties, failureproperties = None, z0=None, z1=None, Sigma = None, Epsilon = None):
        # elasticproperties format: [E1, E2, G12, v12]
        # failureproperties format: [E11f, v21f, msf, R11t, R11c]

        # geometric properties ply
        self.theta = theta
        self.t = t
        self.z0 = None
        self.z1 = None

        # elastic properties
        self.E1 = elasticproperties[0]     # Modulus of elasticity in material direction 1
        self.E2 = elasticproperties[1]     # Modulus of elasticity in material direction 2
        self.G12 = elasticproperties[2]    # Shear modulus in plane 12
        self.v12 = elasticproperties[3]    # Poisson's ratio in plane 12

        #failure properties
        self.E11f = failureproperties[0]   # fiber E11
        self.v21f = failureproperties[1]   # fiber v21
        self.msf = failureproperties[2]    # compensation factor, msf = 1.3 for GFRP, msf = 1.1 for CFRP
        self.R11t = failureproperties[3]   # ultimate tensile stress allowable
        self.R11c = failureproperties[4]   # ultimate compressive stress allowable
        self.Yt = failureproperties[5]
        self.Yc = failureproperties[6]
        self.S = failureproperties[7]
        self.p12 = 0.3

        self.Epsilon = Epsilon
        self.Sigma = Sigma

        self.theta_rad = np.radians(theta)
        m = np.cos(self.theta_rad)  # Cosine of the angle of the material direction 1 with the x-axis
        n = np.sin(self.theta_rad)  # Sine of the angle of the material direction 1 with the x-axis
        # Now let's translate the equations to Python using numpy operations

        self.v21 = self.v12 * self.E2 / self.E1

        Q = 1 - self.v12 * self.v21
        Q11 = self.E1 / Q
        Q22 = self.E2 / Q
        Q12 = self.v12 * self.E2 / Q
        Q66 = self.G12

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
        self.Smatrix = np.linalg.inv(self.Q)

    # We use the following method to carry out stress analysis for the laminate:
    def StressAnalysis(self):
        self.Sigma = np.zeros((3, 1))
        self.Sigma = self.Q @ self.Epsilon # this is the rotated sigma in the xyz frame
        return self.Sigma

    # Finally, we can carry out failure analysis to figure out whether a lamina has failed:
    def FailureAnalysis(self, sigma):
        #first we rotate the stress vector sigma by -theta -> so back into 123 frame
        m = np.cos(-self.theta_rad)
        n = np.sin(-self.theta_rad)
        alfamatrix = np.array([[m, n, 0],
                               [-n, m, 0],
                               [0, 0, 1]])
        sigma123 = alfamatrix @ sigma
        IFFfactor = self.IFF(sigma123)
        FFfactor = self.FF(sigma123)

        if IFFfactor >= 1 or FFfactor >= 1:
            failure = True
        else:
            failure = False
        if IFFfactor >= 1.1 or FFfactor >= 1.1:
            print('failure criteria > 1.1, load increment too big!')

        self.failure = failure
        return failure

    def IFF(self, sigma):
        s1 = sigma[0]
        s2 = sigma[1]
        s6 = sigma[2]

        # We are looking for IFF or inter fiber fracture -> matrix failure due to normal stress in 2 direction
        # or shear stress in 21 or 23 direction (s6). This means we only want to look at this if
        # s2 and s6 are nonzero
        if s2 == 0 or s6 == 0:
            return False
        # Now we have done that first check so we can make the criteria:

        if s2 > 0:
            f = np.sqrt((s6 / self.S) ** 2 + (1 - self.p12 * self.Yt / self.S) ** 2 * (s2 / self.Yt) ** 2) + 0.3 * s2 / self.S
            print('Mode A failure in ply, failure f:', f)

        # Intermediary values:
        s23A = (self.S / (2 * 0.2)) * (np.sqrt(1 + 2 * 0.2 * self.Yc / self.S) - 1)
        p12_minus = 0.25
        p23_minus = p12_minus * s23A / self.S  # p12(-)????
        s12c = self.S * np.sqrt(1 + 2 * p23_minus)

        # figure out mode: B or C
        if abs(s2 / s6) <= s23A / abs(s12c) and abs(s2 / s6) >= 0:
            print('mode B')
            f = (np.sqrt(s6 ** 2 + (p12_minus * s2) ** 2) + p12_minus * s2) / self.S

        else:
            print('mode C')
            f = ((s6/(2*(1 + p23_minus*self.S)))**2 + (s2/self.Yc)**2)*(self.Yc/-s2)
        return f

    def FF(self, sigma):
        s1 = sigma[0]
        s2 = sigma[1]
        if s1 < 0:
            print('lamina is in compression')
            R11 = self.R11t
        else:
            print('lamina is in tension')
            R11 = -1 * self.R11c
        criterion = (1/R11) * (s1 - (self.v21 - self.v21f * self.msf * (self.E1/self.E11f))*s2)
        return criterion

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

    #carry out failure analysis for all lamina in laminate
    def Failurecriteria(self):
        # We need to make sure the lamina have stresses:
        self.StressAnalysis()
        for count, lamina in enumerate(self.laminas):
            lamina.FailureAnalysis(lamina.Sigma)
        pass

#now we test the code:
E1 = 150e9
E2 = 20e9
G12 = 5e9
v12 = 0.3
elasticproperties = [E1, E2, G12, v12]

E11f = 500e9
v21f = 0
msf = 1.1
R11t = 1000e6
R11c = 800e6
yt = 100e6
yc = 100e6
S = 100e6
failureproperties = [E11f, v21f, msf, R11t, R11c, yt, yc, S]
s0 = Lamina(0.0002, 0, elasticproperties, failureproperties)
s1 = Lamina(0.0002, 0, elasticproperties, failureproperties)
s2 = Lamina(0.0002, 0, elasticproperties, failureproperties)
s3 = Lamina(0.0002, 0, elasticproperties, failureproperties)

# create the laminas list, for the laminate function:
laminas = [s0, s1, s2, s3]

# creating the laminate object:
laminate = Laminate(laminas)

# now we can apply loads to the laminate: (in Newton per meter or newtonmeter per meter)
laminate.Loads = np.array([[10000], [0], [0], [0], [0], [0]])

laminate.Failurecriteria()

