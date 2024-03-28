import numpy as np
class Lamina:
    def __init__(self, t, theta, elasticproperties, statisticalproperties = None, failureproperties = None, z0=None, z1=None, Sigma = None, Epsilon = None, FailureState = 0):
        # elasticproperties format:     [E1, E2, G12, v12]
        # failureproperties format:     [E11f, v21f, msf, R11t, R11c]
        # statisticalproperties format: [E1, E2, v12, G12, Xt, Xc, Yt, Yc, S] all standard deviations
        # Damagrpropagation:
        # [0.1, 0, 0.1, 0.1] for failurestate 1
        # [0, 0, 0, 0] for failurestate 2
        self.DamageProgression = np.array([[1, 1, 1, 1],
                                           [0.1, 1e-10, 0.1, 0.1],
                                           [1e-10, 1e-10, 1e-10, 1e-10]])

        # geometric properties ply
        self.theta = theta
        self.t = t
        self.z0 = None
        self.z1 = None

        # We save the original elastic properties:
        self.ElasticProperties = elasticproperties

        # Then we save the current elastic properties based on the failure state:
        self.FailureState = FailureState
        self.ElasticPropertiesCurrent = self.DamageProgression[self.FailureState] * self.ElasticProperties

        # Now we define the elastic properties in case we want to call them right after initialisation
        self.elasticproperties = elasticproperties
        self.E1 = self.ElasticPropertiesCurrent[0]     # Modulus of elasticity in material direction 1
        self.E2 = self.ElasticPropertiesCurrent[1]     # Modulus of elasticity in material direction 2
        self.G12 = self.ElasticPropertiesCurrent[2]    # Shear modulus in plane 12
        self.v12 = self.ElasticPropertiesCurrent[3]    # Poisson's ratio in plane 12

        # failure properties
        self.failureproperties = failureproperties
        self.E11f = failureproperties[0]   # fiber E11
        self.v21f = failureproperties[1]   # fiber v21
        self.msf = failureproperties[2]    # compensation factor, msf = 1.3 for GFRP, msf = 1.1 for CFRP
        self.R11t = failureproperties[3]   # ultimate tensile stress allowable
        self.R11c = failureproperties[4]   # ultimate compressive stress allowable
        self.Yt = failureproperties[5]
        self.Yc = failureproperties[6]
        self.S = failureproperties[7]
        self.p12 = 0.3

        # Failure states:
        self.Epsilon = Epsilon
        self.Sigma = Sigma
        self.FailureStresses = []

        self.theta_rad = np.radians(theta)
        self.CalculateQS()

    def CalculateQS(self):
        # # If the failure state is 1 or above for now just remove the elastic properties

        # Based on the failurestate, we should take different values for the elastic properties:
        self.ElasticPropertiesCurrent = self.DamageProgression[self.FailureState] * self.ElasticProperties

        # Now we update the elastic property attributes for later use in failure criteria
        self.E1 = self.ElasticPropertiesCurrent[0]
        self.E2 = self.ElasticPropertiesCurrent[1]
        self.G12 = self.ElasticPropertiesCurrent[2]
        self.v12 = self.ElasticPropertiesCurrent[3]
        self.v21 = self.v12 * self.E2 / self.E1

        m = np.cos(self.theta_rad)  # Cosine of the angle of the material direction 1 with the x-axis
        n = np.sin(self.theta_rad)  # Sine of the angle of the material direction 1 with the x-axis

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
        try:
            self.Smatrix = np.linalg.inv(self.Q)
        except:
            print("Singular Q matrix, does this mean it's zero?")

    # We use the following method to carry out stress analysis for the laminate:
    def StressAnalysis(self):
        self.Sigma = np.zeros((3, 1))
        self.Sigma = self.Q @ self.Epsilon # this is the rotated sigma in the xyz frame
        return self.Sigma

    # Finally, we can carry out failure analysis to figure out whether a lamina has failed:
    def FailureAnalysis(self, sigma):
        # first we rotate the stress vector sigma by -theta -> so back into 123 frame
        m = np.cos(-self.theta_rad)
        n = np.sin(-self.theta_rad)
        alfamatrix = np.array([[m**2, n**2, 2*m*n],
                               [n**2, m**2, -2*m*n],
                               [-m*n, m*n, m**2 - n**2]])
        sigma123 = alfamatrix @ sigma
        IFFfactor = self.IFF(sigma123)
        FFfactor = self.FF(sigma123)

        if IFFfactor >= 1 or FFfactor >= 1:
            failure = 1

            # We want to save the stesses at failure of a ply:
            self.FailureStresses.append(sigma123)
        else:
            failure = 0

        # Modify the failure state based on whether failure is achieved or not
        self.FailureState += failure

        # We'll return the actual factor as well, so we can use it to find the next loadstep:
        return failure, IFFfactor, FFfactor

    # This function does inter fiber failure analysis
    # It does not alter anything about the failure state of the material
    def IFF(self, sigma):
        s2 = sigma[1]
        s6 = sigma[2]

        # Intermediary values:
        p12_minus = 0.25
        s23A = (self.S / (2 * p12_minus)) * (np.sqrt(1 + 2 * p12_minus * self.Yc / self.S) - 1)
        p23_minus = p12_minus * s23A / self.S  # p12(-)????
        s12c = self.S * np.sqrt(1 + 2 * p23_minus)

        # Mode A:
        if s2 >= 0:
            f = np.sqrt((s6 / self.S)**2 + (1 - self.p12 * self.Yt / self.S)**2 * (s2 / self.Yt)**2) + 0.3 * s2 / self.S
        # Now if s2 < 0 then it could be either mode b or mode c, which we check:
        elif abs(s2 / (abs(s6) + 1e-11)) <= (s23A / abs(s12c)) and abs(s2 / (abs(s6) + 1e-11)) >= 0:
            # In the above line the following code: (abs(s6) + 1e-11) functions to prevent a div by zero error
            f = (np.sqrt(s6 ** 2 + (p12_minus * s2) ** 2) + p12_minus * s2) / self.S
        else:
            term1 = (s6/(2*(1 + p23_minus)* self.S))
            term2 = s2/self.Yc
            f = (term1**2 + term2**2)*(self.Yc/-s2)
        return f


    # This function does fiber failure analysis
    # It does not alter anything about the failure state of the material
    def FF(self, sigma):
        s1 = sigma[0]
        s2 = sigma[1]
        if s1 < 0:
            R11 = -self.R11c
        else:
            R11 = self.R11t
        f = (1 / R11) * (s1 - (self.v21 - self.v21f * self.msf * (self.E1 / self.E11f)) * s2)
        return f

