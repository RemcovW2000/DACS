import numpy as np
import scipy.optimize as opt

from structural_entity import StructuralEntity
class Lamina(StructuralEntity):
    """
    Represents a single ply of a composite laminate.

    This object is used in the laminate class, which forms a stack of Lamina objects ordered in the 'laminas' list. The
    lamina object can obtain load intensities (N/mm) either from a user directly (rare) or from a Laminate object.

    Attributes:
        t (float): Ply thickness in mm.
        theta_ (float): Ply orientation angle in degrees, where 0 describes the lamina to be aligned with the x-axis.
        elastic_properties (list of float): Elastic constants [E1, E2, G12, v12].
        failure_state (int): 0 for no failure, 1 for inter-fiber failure, 2 for fiber failure.
        elastic_properties_current (list of float): Updated elastic properties after failure.

    Methods:
        CalculateQS(): Calculates stiffness (Q) and compliance (S) matrices.
        StressAnalysis(): Computes stress state based on current strain state.
        FailureAnalysis(sigma): Determines the failure state of the ply.
    """
    def __init__(self, t, theta, elastic_properties, failureproperties = None, rho = None, z0=None, z1=None, Sigma = None, Epsilon = None, failure_state = 0):
        super().__init__('lamina')
        # elastic_properties format:     [E1, E2, G12, v12]
        # failureproperties format:     [E11f, v21f, msf, R11t, R11c]
        # statisticalproperties format: [E1, E2, v12, G12, Xt, Xc, Yt, Yc, S] all standard deviations
        # Damagepropagation:
        # [0, 0.1, 0.1, 0.1] for failure_state 1
        # [0, 0, 0, 0] for failure_state 2
        # The failure state of the lamina is 0 for an undamaged ply, and 1 for a ply that has sustained-
        # an IFF. If the ply sustains a FF or more than one IFF, the failure state is 2.
        self.damage_progression_ = np.array([[1, 1, 1, 1],
                                           [1, 0.1, 0.1, 0.1],
                                           [1e-10, 1e-10, 1e-10, 1e-10]])

        # The following line determines whether the max stress criteria is used or puck:
        self.max_stress = False

        # geometric properties ply
        self.theta_ = theta
        self.t = t

        # Z coordinate of plys is assigned by the laminate class
        self.z0 = None
        self.z1 = None

        # Now we define the elastic properties in case we want to call them right after initialisation
        self.elastic_properties = elastic_properties
        self.elastic_properties_current = elastic_properties
        self.E1 = self.elastic_properties_current[0]     # Modulus of elasticity in material direction 1
        self.E2 = self.elastic_properties_current[1]     # Modulus of elasticity in material direction 2
        self.G12 = self.elastic_properties_current[2]    # Shear modulus in plane 12
        self.v12 = self.elastic_properties_current[3]    # Poisson's ratio in plane 12

        self.rho = rho

        # Then we save the current elastic properties based on the failure state:
        self.failure_state = failure_state
        self.elastic_properties_current = self.damage_progression_[self.failure_state] * self.elastic_properties

        # failure properties
        self.failureproperties = failureproperties
        self.E11f = failureproperties[0]    # fiber E11
        self.v21f = failureproperties[1]    # fiber v21
        self.msf = failureproperties[2]     # compensation factor, msf = 1.3 for GFRP, msf = 1.1 for CFRP
        self.R11t = failureproperties[3]    # ultimate tensile stress allowable
        self.R11c = failureproperties[4]    # ultimate compressive stress allowable
        self.Yt = failureproperties[5]      # Matrix tensile strength
        self.Yc = failureproperties[6]      # Matrix compressive strength
        self.S = failureproperties[7]       # Matrix shear strength
        self.p12 = 0.3                      # p12 for CFRP is 0.3

        # Lamina states:
        self.Epsilon = Epsilon              # Strain state
        self.Sigma = Sigma                  # Stress state
        self.FailureStresses = []           # List containing stresses at failure

        # Hard coded for now
        # TODO put this into MP file
        self.Sxz = 92 # Shear strength in xz (in direction of fibers)
        self.Syz = 73 # Shear strength in yz (transversely to fibers)


        # Lastly we calculate the Q and S (stifness and compliance) matrices in the initialisation
        self.calculate_QS()

    def calculate_QS(self):
        # Based on the failure_state, we should take different values for the elastic properties:
        if self.failure_state <= 2:
            self.elastic_properties_current = self.damage_progression_[self.failure_state] * self.elastic_properties
        else:
            # If failure state is 2 or greater we automatically overrule the damage_progression_ matrix-
            # and apply zero elastic properties
            self.elastic_properties_current = [1e-10, 1e-10, 1e-10, 1e-10]

        # Now we update the elastic property attributes for later use in failure criteria
        self.E1 = self.elastic_properties_current[0]
        self.E2 = self.elastic_properties_current[1]
        self.G12 = self.elastic_properties_current[2]
        self.v12 = self.elastic_properties_current[3]
        self.v21 = self.v12 * self.E2 / self.E1

        # Now calculating values for the Q matrix:
        m = np.cos(np.deg2rad(self.theta_))  # Cosine of the angle of the material direction 1 with the x-axis
        n = np.sin(np.deg2rad(self.theta_))  # Sine of the angle of the material direction 1 with the x-axis

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

        # Generating Q matrix
        self.Q = np.array([[Qxx, Qxy, Qxs],
                           [Qxy, Qyy, Qys],
                           [Qxs, Qys, Qss]])
        # This Q is the ROTATED (global) q matrix, we can invert it to get the compliance matrix

        try:
            self.Smatrix = np.linalg.inv(self.Q)
        except:
            # it may be that the Q matrix is singular, or inversion does not work so we use the except statement
            print("Singular Q matrix")

    # We use the following method to carry out stress analysis for the laminate:
    def stress_analysis(self):
        self.Sigma = np.zeros((3, 1))
        self.Sigma = self.Q @ self.Epsilon # this is the rotated sigma in the xyz frame
        return self.Sigma

    # Finally, we can carry out failure analysis to figure out whether a lamina has failed:
    def failure_analysis(self):
        super().failure_analysis()
        sigma = self.Sigma
        # first we rotate the stress vector sigma by theta_ -> so back into 123 frame
        m = np.cos(np.deg2rad(self.theta_))
        n = np.sin(np.deg2rad(self.theta_))
        alfamatrix = np.array([[m**2, n**2, 2*m*n],
                               [n**2, m**2, -2*m*n],
                               [-m*n, m*n, m**2 - n**2]])
        sigma123 = alfamatrix @ sigma

        # Choose whether to use max stress or puck:
        if self.max_stress == True:
            IFF_factor = self.IFFmax(sigma123)
            FF_factor = self.FFmax(sigma123)
        else:
            IFF_factor = self.IFF(sigma123)
            FF_factor = self.FF(sigma123)

        self.set_failure_indicator('inter_fiber_failure', IFF_factor)
        self.set_failure_indicator('fiber_failure', FF_factor)

        # Now we check whether either of the criteria factors are greater than one
        if IFF_factor >= 1:
            failure = 1
            self.FailureStresses.append(sigma123)
        elif FF_factor >= 1:
            failure = 2
            # We want to save the stesses at failure of a ply:
            self.FailureStresses.append(sigma123)
        else:
            failure = 0

        # Modify the failure state based on whether failure is achieved or not
        self.failure_state += failure

        # We'll return the actual factor as well, so we can use it to calculate the next loadstep:
        return failure, IFF_factor, FF_factor

    # This function does inter fiber failure analysis
    # It does not alter anything about the failure state of the material, simpy returns the factor
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
    # It does not alter anything about the failure state of the material, simply returns the factor
    def FF(self, sigma):
        s1 = sigma[0]
        s2 = sigma[1]

        # Choosing which strength to use, compression or tension:
        if s1 < 0:
            R11 = -self.R11c
        else:
            R11 = self.R11t

        f = (1 / R11) * (s1 - (self.v21 - self.v21f * self.msf * (self.E1 / self.E11f)) * s2)
        return f

    # This is for the max stress:
    def FFmax(self, sigma):
        s1 = sigma[0]
        if s1 < 0:
            R11 = -self.R11c
        else:
            R11 = self.R11t
        f = np.abs(s1)/np.abs(R11)
        return f

    # Again max stress
    def IFFmax(self, sigma):
        s2 = sigma[1]
        s6 = sigma[2]
        if s2 < 0:
            y = -self.Yc
        else:
            y = self.Yt

        ftrans = np.abs(s2)/np.abs(y)
        fshear = np.abs(s6)/np.abs(self.S)
        return max(ftrans, fshear)

    # This code was used for question 1b: only for verifying the stresses found in the graphs generated with code v1
    def calculate_principal_sigma_epsilon(self):
        # The stress analysis is a prerequisite for this as we need to know the stresses
        self.stress_analysis()

        #Now we make a rotation matrix to rotate the stresses into the principal frame
        m = np.cos(np.deg2rad(self.theta_))
        n = np.sin(np.deg2rad(self.theta_))

        # Multiply global stresses by rotation matrix:
        AlfaStress = np.array([[m**2, n**2, 2*m*n],
                               [n**2, m**2, -2*m*n],
                               [-m*n, m*n, m**2 - n**2]])
        Sigma123 = AlfaStress @ self.Sigma

        # Now do the same for strains
        AlfaStrain = np.array([[m**2, n**2, m*n],
                               [n**2, m**2, -m*n],
                               [-2*m*n, 2*m*n, m**2 - n**2]])
        Epsilon123 = AlfaStrain @ self.Epsilon

        # Return the principal stresses and strains:
        return Sigma123, Epsilon123

    def delamination_analysis(self, Taurz0, Taurz1, azimuth):
        # Taurz0 = shear stress at bottom of ply
        # Taurz1 = shear stress at top of ply

        bottomdelamination = False
        topdelamination = False

        # angle between direction of stress and the 123 axis system:
        psi = self.theta_-azimuth

        Tau1z0 = Taurz0 * np.cos(np.deg2rad(psi))
        Tau2z0 = Taurz0 * np.sin(np.deg2rad(psi))

        Tau1z1 = Taurz1 * np.cos(np.deg2rad(psi))
        Tau2z1 = Taurz1 * np.sin(np.deg2rad(psi))

        # Now we can use max stress:
        if abs(Tau1z0/self.Sxz) >= 1 or abs(Tau2z0/self.Syz) >= 1:
            bottomdelamination = True

        if abs(Tau1z1/self.Sxz) >= 1 or abs(Tau2z1/self.Syz) >= 1:
            topdelamination = True
        # else:
            # print('stresses on top surface: ', Tau1z1, Tau2z1)

        return bottomdelamination, topdelamination

    def Taucrit(self, azimuth):
        """
        This function calculates the interlaminar shear strength in a given direction
        :param azimuth:
        :returns: scalar value for interlaminar shear strength
        """
        psi = self.theta_ - azimuth

        def Find0(Taur):
            term1 = Taur*np.cos(np.deg2rad(psi))/self.Sxz
            term2 = Taur*np.sin(np.deg2rad(psi))/self.Syz

            Fi = term1**2 + term2**2
            return Fi-1

        # define bounds for the solution of Taur -> cannot be outside the two values for S
        a = self.Syz
        b = self.Sxz

        # Find the root within the given bounds
        Taucrit = opt.brentq(Find0, a, b)

        return Taucrit

    def calculate_weight_per_A(self):
        return self.t*self.rho