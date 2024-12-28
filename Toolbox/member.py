# External packages
from tqdm import tqdm
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp

# Local imports
from Toolbox.damaged_region import *
from Toolbox.curved_plate_data import k_function
from Toolbox.laminate import laminate_builder
from structural_entity import StructuralEntity

class Member(StructuralEntity):
    def __init__(self, panel, loads = [0, 0, 0, 0, 0, 0], a = 300, b = 200):
        """"
        This class makes a member: either a sandwich - or laminate object, with certain dimensions:
        1. Width -> b
        In direction of y-axis
        2. Depth -> a
        in direction of x-axis
        3. thickness -> h
        determined by layup

        Damage analysis:
        1. run self.Fimpact() -> calculates force at impact
        """
        super().__init__('member')
        self.panel = panel              # could be either laminate or sandwich, both should work
        self.Loads = loads              # total loads, not intensity!!!
        self.submember = None           # this is a submember for adding local reinforcement.
        self.submember_start = None     # x coordinate of the start of the submember, in airfoil FOR
        self.submember_end = None       # x coordinate of end of the submember, in airfoil FOR

        # this submember is just an instance of the member class

        self.startcoord = [] # from leading to trailing edge for wing
        self.endcoord = [] # from leading to trailing edge for wing
        self.R = None
        if not self.startcoord:
            self.b = b # panel width
        self.a = a  # panel depth
        self.curvealongx = None
        self.booms = None
        self.segments = None
        self.arclength = None

        self.xbar = None
        self.ybar = None

        self.weight_per_b_ = None
        # ----------------------------------------------------------------
        # Assign failure indicators as None:
        self.BucklingFI = None
        self.panel_FI = None

        # ----------------------------------------------------------------
        # Dacs II: impact
        self.R_impactor = 6.35
        self.E_impactor = 200000  # MPa
        self.v_impactor = 0.28

        self.BVID_energy = 0.113*2000/25.4 # Joules per inch -> mm thickness
        # ----------------------------------------------------------------

    @property
    def child_objects(self):
        return [self.panel, self.submember]

    def get_Ex(self, x):
        '''
        Function obtains the youngs modulus in the x direction for a member.

        This function was implemented to correctly obtain the Ex for any member, even if it has a submember. As
        submembers are defined by two extreme x values within which the submember is located, it takes an x coordinate
        as an input.
        :param x:
        :return: Ex:
        '''
        # function should return Ex based on value for x, whether it's in the submember or not.
        # reference point for x = 0 is the global reference point
        if self.submember:
            if x >= self.submember_start and x <= self.submember_end:
                Ex = self.submember.panel.Ex
            else:
                Ex = self.panel.Ex
        else:
            Ex = self.panel.Ex
        return Ex

    def get_h(self, x):
        '''
        Function obtains the thickness in the x direction for a member.

        This function was implemented to correctly obtain the thickness h for any member, even if it has a submember. As
        submembers are defined by two extreme x values within which the submember is located, it takes an x coordinate
        as an input.
        :param x:
        :return: h:
        '''
        # function should return h based on value for x, whether it's in the submember or not.
        # submember is a member type object!
        # reference point for x = 0 is the global reference point
        if self.submember:
            if x >= self.submember_start and x <= self.submember_end:
                h = self.submember.panel.h
            else:
                h = self.panel.h
        else:
            h = self.panel.h
        return h

    def failure_analysis(self):
        '''
        This function carries out a failure analysis of the member specifically for a cross-sectional analysis (CSA).

        The function uses the present booms in order to obtain the stresses of the laminate, and assigns these stresses
        to the correct laminate (submember or full member)

        most conservative is to take the max shear stress and max normal stress in each member as stress state for buckling
        as well as FPF

        more accurate is to check FPF at every segment/boom stress combination
        '''
        super().failure_analysis()
        # -------------------------------------------------------------------------
        # Calculating load cases:
        # -------------------------------------------------------------------------
        if self.booms:
            if self.submember:
                submember_booms = [boom for boom in self.booms if self.submember_start < boom.location[0] + self.xbar < self.submember_end]
                booms =[boom for boom in self.booms if self.submember_start > boom.location[0] + self.xbar or self.submember_end < boom.location[0] + self.xbar]
                # A segment is considered in the submember if one of the booms is in the segment, such that the potentially
                # higher stress is not counted into the main member

                submember_segments = [segment for segment in self.segments if self.submember_start < min(segment.p1[0]+self.xbar, segment.p2[0]+self.xbar) and self.submember_end > max(segment.p1[0]+self.xbar, segment.p2[0]+self.xbar)]
                segments = [segment for segment in self.segments if self.submember_start > max(segment.p1[0] + self.xbar, segment.p2[0]+self.xbar) or self.submember_end < min(segment.p1[0]+self.xbar,segment.p2[0]+self.xbar)]

                self.submember.booms = submember_booms
                self.submember.segments = submember_segments

                submember_max_FI = self.submember.failure_analysis()
            else:
                booms = self.booms
                segments = self.segments

            normal_force_intensities = [(boom.Sigmax * self.get_h(boom.location[0])) for boom in booms]

            # we must find the most critical stress state, which could in fact be multiple! positive or negative.
            max_normal_force_intensity = max(normal_force_intensities)
            min_normal_force_intensity = min(normal_force_intensities)

            shear_flows = [segment.qs for segment in segments] # N/mm

            # again for shear:
            max_shear_flow = max(shear_flows)
            min_shear_flow = min(shear_flows)

            # make several load cases: combining min with max etc.:
            extreme_load_case_1 = [max_normal_force_intensity, 0, max_shear_flow, 0, 0 ,0]
            extreme_load_case_2 = [min_normal_force_intensity, 0, min_shear_flow, 0, 0 ,0]
            extreme_load_case_3 = [max_normal_force_intensity, 0, min_shear_flow, 0, 0 ,0]
            extreme_load_case_4 = [min_normal_force_intensity, 0, max_shear_flow, 0, 0 ,0]
            loadcases = [extreme_load_case_1, extreme_load_case_2, extreme_load_case_3, extreme_load_case_4]

            # -------------------------------------------------------------------------
            # Applying load cases to failure modes:
            # -------------------------------------------------------------------------

            # first ply failure:
            failureindicators_panel = []
            for loadcase in loadcases:
                self.panel.Loads = loadcase
                failureindicators_panel.append(self.panel.failure_analysis())

            self.panel_FI = max(failureindicators_panel)

            # member buckling:
            index_max_abs_shear = max(enumerate(shear_flows), key=lambda x: abs(x[1]))[0]

            max_normal = min(normal_force_intensities)
            max_abs_shear = shear_flows[index_max_abs_shear]

            self.Loads = [max_normal * self.b, 0, max_abs_shear * self.a, 0, 0 ,0]

            # Curvature must have been assigned in airfoil class!
            Buckling_FI = self.buckling_analysis()
            self.BucklingFI = Buckling_FI
            FI_list = failureindicators_panel + [Buckling_FI]
            print(failureindicators_panel, Buckling_FI)
            max_FI = max(FI_list)
            if self.submember:
                FI_list.append(submember_max_FI)
                max_FI = max(FI_list)

            self.set_failure_indicator('buckling', Buckling_FI)
            self.set_failure_indicator('child', max([max(value for key, value in
                                                         child_object.failure_indicators.items() if
                                                         isinstance(value, (int, float))) for child_object in
                                                     self.child_objects if child_object]))
        else:
            raise TypeError('No booms exist, code for use of members outside of cross sectional analysis must be'
                            ' finised, add a way to run both buckling and first ply failure analysis using assigned'
                            'loads instead of loadcases derived from stress state in booms.')
        return max(value for key, value in self.failure_indicators.items() if isinstance(value, (int, float)))

    def buckling_analysis(self):
        '''
        Base function to perform the buckling analysis for the member. Loads must already be assigned!

        It should work even if curvature is zero, or certain parameters are not assigned such as R (radius)

        :return:
        '''
        # find the right scaling factor:
        if self.R == None:
            print('R not assigned!')
            # R not assigned -> use flat plate eq:
            scalingfactor = 1
        elif self.b**2/self.R < 1:
            # too flat for curved equations: use flat plate eq:
            scalingfactor = 1
        elif self.b**2/self.R > 100:
            # use equations for full cylinder!
            print('full cylinder equations needed, not programmed! Outcome is false, (and is assigned: Ncrit = 0)')
            scalingfactor = 1
            pass
        else:
            # use curved plate equation:
            Zb = np.sqrt(1-self.panel.vxy*self.panel.vyx)*self.b**2/self.R
            Rovert = self.R/self.panel.h
            scalingfactor = k_function(Zb, Rovert)

        # ------------------------------------------------------------------------------------
        # find the failure indicator:
        if self.Loads[0] > 0:
            # TODO: add a more nuanced way to check this! shear load could still cause buckling!
            FI = 0
        else:
            Ncritnormal = self.combined_load_buckling_Ncrit(self.Loads)
            Ncrit = scalingfactor * Ncritnormal
            FI = self.Loads[0]/Ncrit
        return FI

    def combined_load_buckling_Ncrit(self, Loads):
        D11 = self.panel.ABD_matrix[3, 3]
        D12 = self.panel.ABD_matrix[3, 4]
        D22 = self.panel.ABD_matrix[4, 4]
        D66 = self.panel.ABD_matrix[5, 5]

        # we check the formula for many m (# of half waves in this buckling mode)
        Ncritnorm = (np.pi / self.a) ** 2 * (
                    D11 + 2 * (D12 + 2 * D66) * (self.a / self.b) ** 2 + D22 * (
                        self.a / self.b) ** 4)

        Nx = Loads[0] / self.b
        Nxy = Loads[2] / self.b
        k = Nxy/Nx
        denominator = 2 - (8192 * self.a**2 * k**2)/(81 * self.b**2 * np.pi**4)

        term2 = 5 - np.sqrt(9 + (65536*self.a**2*k**2)/(81*np.pi**4*self.b**2))
        Ncrit = abs((Ncritnorm/denominator)*term2)

        # Now obtain the scaling factor from the panel:
        Ncrit = self.panel.buckling_scaling_factor(Ncrit)
        return Ncrit

    def calculate_b(self):
        x1, y1 = self.startcoord
        x2, y2 = self.endcoord
        self.b = np.sqrt((x2-x1)**2 + (y2-y1)**2)
        return self.b

    def calculate_Zb(self):
        if self.R == None:
            self.Zb = 0
        else:
            self.Zb = np.sqrt(1-self.panel.vxy*self.panel.vyx)*self.b**2/(self.R*self.panel.h)
        return self.Zb

    def shear_load_buckling_FI(self):
        D11 = self.panel.ABD_matrix[3, 3]
        D12 = self.panel.ABD_matrix[3, 4]
        D22 = self.panel.ABD_matrix[4, 4]
        D66 = self.panel.ABD_matrix[5, 5]

        term1 = (9 * np.pi ** 4 * self.b) / (32 * self.a ** 3)
        term2 = (D11 + 2 * (D12 + 2 * D66)) * (self.a ** 2 / self.b ** 2)
        term3 = D22 * (self.a ** 4 / self.b ** 4)
        Nxypanel = 0.78 * term1 * (term2 + term3)
        return Nxypanel

    def calculate_panel_FI(self):
        """
        returns critical load in loading ratio of applied load. Load must be divided by width!
        :return:
        """
        # Divide the x load by the width -> b
        a = self.a
        b = self.b
        Loads = self.Loads
        Nx = Loads[0]/b
        Ny = Loads[1]/a
        Ns = Loads[2]/b
        Mx = Loads[3]/b
        My = Loads[4]/a
        Mz = Loads[5]/b # not sure abt this one!

        # Now we have load intensities:
        Loadintensities = [Nx, Ny, Ns, Mx, My, Mz]
        self.panel.Loads = Loadintensities
        Ncrit = self.panel.Ncrit()
        Fx = Ncrit[0]*b
        Fy = Ncrit[1]*a
        Fs = Ncrit[2]*b
        Mx = Ncrit[3]*b
        My = Ncrit[4]*a
        Mz = Ncrit[5]*b

        Fcrit = [Fx, Fy, Fs, Mx, My, Mz]
        return Fcrit

    def normal_load_buckling_FI(self):
        """
        In the current loading ratio, what is the critical load that would cause buckling?
        :return:
        """
        D11 = self.panel.ABD_matrix[3, 3]
        D12 = self.panel.ABD_matrix[3, 4]
        D22 = self.panel.ABD_matrix[4, 4]
        D66 = self.panel.ABD_matrix[5, 5]

        Fcrits = []
        for m in range(1, 6):
            Fcrit = self.b * ((np.pi / self.a) ** 2 * (D11*m**2 + 2 * (D12 + 2 * D66) * (self.a/self.b) ** 2 + (D22/m**2) * (self.a/self.b) ** 4))
            Fcrits.append(Fcrit)

        Fcrit_min = min(Fcrits)
        min_index = Fcrits.index(Fcrit_min)
        min_mode = min_index+1
        return Fcrit_min

    def calculate_arc_length(self):
        '''
        Calculates the arc length of the member in the xz plane in wing level analysis.

        :return: float
        '''
        return np.sum([np.sqrt((segment.p2[0]-segment.p1[0])**2 + (segment.p2[1]-segment.p1[1])**2) for segment in self.segments])

    def calculate_weight_per_b(self):
        '''
        Calculates the weight per unit span of the member given the arc length of the member and the weight of the panel

        :return: None
        '''

        if self.submember:
            sub_member_segments = [segment for segment in self.segments if self.submember_start < min(segment.p1[0]+self.xbar, segment.p2[0]+self.xbar) and self.submember_end > max(segment.p1[0]+self.xbar, segment.p2[0]+self.xbar)]
            segments = [segment for segment in self.segments if self.submember_start > max(segment.p1[0] + self.xbar, segment.p2[0]+self.xbar) or self.submember_end < min(segment.p1[0]+self.xbar,segment.p2[0]+self.xbar)]

            normal_member_arc_length = np.sum([np.sqrt((segment.p2[0] - segment.p1[0]) ** 2 + (segment.p2[1] - segment.p1[1]) ** 2) for segment in
                    segments])
            normal_panel_weight_per_A = self.panel.calculate_weight_per_A() * normal_member_arc_length

            sub_member_arc_length = np.sum([np.sqrt((segment.p2[0] - segment.p1[0]) ** 2 + (segment.p2[1] - segment.p1[1]) ** 2) for segment in
                    sub_member_segments])
            sub_member_panel_weight_per_A = self.submember.panel.calculate_weight_per_A()*sub_member_arc_length
            weight_per_b = normal_panel_weight_per_A + sub_member_panel_weight_per_A

        else:
            normal_member_arc_length = np.sum([np.sqrt((segment.p2[0] - segment.p1[0]) ** 2 + (segment.p2[1] - segment.p1[1]) ** 2) for segment in
                    self.segments])
            normal_panel_weight_per_A = self.panel.calculate_weight_per_A() * normal_member_arc_length
            weight_per_b = normal_panel_weight_per_A

        self.weight_per_b_ = weight_per_b
        return weight_per_b

    def deflection_single_term(self, F, m, n, xo, yo, x, y):
        a = self.a
        b = self.b

        ABD = self.panel.ABD_matrix

        D11 = ABD[3, 3]
        D12 = ABD[3, 4]
        D66 = ABD[5, 5]
        D22 = ABD[4, 4]

        # Calculate components of the numerator
        term1 = np.sin(m * np.pi * xo / a)
        term2 = np.sin(n * np.pi * yo / b)
        term3 = np.sin(m * np.pi * x / a)
        term4 = np.sin(n * np.pi * y / b)
        numerator = 4 * (F / (a * b)) * term1 * term2 * term3 * term4

        # Calculate components of the denominator
        term5 = D11 * (m * np.pi / a) ** 4
        term6 = 2 * (D12 + 2 * D66) * (m ** 2 * n ** 2 * np.pi ** 4 / (a ** 2 * b ** 2))
        term7 = D22 * (n * np.pi / b) ** 4
        denominator = term5 + term6 + term7

        return numerator / denominator

    def compute_deflection(self, F, xo, yo, x, y):
        # Number of terms in fourier series hardcoded for now:
        max_m = 25
        max_n = 25

        result = 0
        for m in range(1, max_m + 1):
            for n in range(1, max_n + 1):
                result += self.deflection_single_term(F, m, n, xo, yo, x, y)
        return result


# ----------------------------------------------------------------------------------------------------------------------
#       DACS II impact damage:
#       Find the force:
# ----------------------------------------------------------------------------------------------------------------------
    def E_deflection(self, F, xo, yo):
        # For the deflection energy we analyse the work done at location of impact so
        # we need the deflection and force at place of impact:
        x = xo
        y = yo

        wmax = self.compute_deflection(F, xo, yo, x, y)
        E_deflection = F*wmax/2
        return E_deflection

    def E_indentation(self, F):
        k = self.k_stiffness(self.R_impactor)
        deltamax = self.calculate_deltamax(F, k)
        E_indentation =0.4*k*(deltamax)**2.5
        return E_indentation

    def calculate_deltamax(self, F, k):
        return (F/k)**(2/3)

    def calculate_Rc(self):
        Rc = np.sqrt(self.R_impactor**2 - (self.R_impactor-self.deltamax)**2)
        return Rc

    def k_stiffness(self, R):
        # Properties of impactor:
        vi = self.v_impactor
        Ei = self.E_impactor

        K1 = (1-vi**2)/Ei

        # ABD = self.panel.ABD_matrix
        # A11 = ABD[0,0]
        # A22 = ABD[1,1]
        # A12 = ABD[0,1]
        #
        # # Ask how to calculate this:
        # Gzr = 4500 # PLACEHOLDER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #
        # term1 = np.sqrt(A22)
        # inner_term1 = np.sqrt(A11 * A22) + Gzr
        # inner_term2 = A12 + Gzr
        # numerator = term1 * np.sqrt(inner_term1 ** 2 - inner_term2 ** 2)
        #
        # # Calculate components of the denominator
        # denominator = 2 * np.pi * np.sqrt(Gzr) * (A11 * A22 - A12 ** 2)
        # K2 = numerator / denominator
        vo = (0.3+0.45)/2
        K2 = (1-vo**2)/11.2e3
        k = 4*np.sqrt(R)/(3*np.pi*(K1 + K2))
        return k

    def E_impact_estimate(self, F, xo, yo):
        return self.E_indentation(F) + self.E_deflection(F, xo, yo)

    def calculate_impactForce(self, xo, yo, tol=1e-4, max_iter=1000):
        # This function sets the attribute Fimpact and also returns the force
        Eimpact = self.BVID_energy*self.panel.h*1000
        print('BVID energy is:', np.round(Eimpact/1000, 1), 'Joules')

        # Let's also save the impact energy for which
        self.Eimpact = Eimpact

        F = 1000 # Initial guess for F
        for i in range(max_iter):
            E_est = self.E_impact_estimate(F, xo, yo)
            diff = E_est - Eimpact

            if abs(diff) < tol:
                # Set impact force attribute for later use:
                self.Fimpact = F

                # Set delta max attribute for later use:
                self.deltamax = self.calculate_deltamax(F, self.k_stiffness(self.R_impactor))
                return F  # Found the solution within tolerance

            # Newton-Raphson update step
            dE_dF = (self.E_impact_estimate(F + tol, xo, yo) - E_est) / tol  # Numerical derivative
            F = F - diff / dE_dF

            if F < 0:
                F = tol  # Ensure F stays positive, you can choose a small positive number

        raise ValueError("Failed to converge to a solution")

    def Taurz(self, r, z):
        # First run the following:
        # 2. self.calculate_impactForce()
        # We don't run this everytime because it's time consuming!
        Ftotal = self.Fimpact
        Rc = self.calculate_Rc()

        if r < Rc:
            if r == 0:
                r = 1e-10
            a = 1/Rc**2
            term1 = Ftotal
            term2 = (1-a*r**2)**1.5
            Fr = term1*(1 - term2)
            Tauavg = Fr/(2*np.pi*r*self.h)
        elif r >= Rc:
            Tauavg = Ftotal/(2*np.pi*r*self.h)

        Taumax = 1.5*Tauavg

        # Now use Tauavg to find the
        Taurz = Taumax*(1-((2*z)**2)/self.h**2)
        return Taurz

    def delamination_analysis(self, azimuth, rmax):
        """
        Checks delaminations at a specific angle azimuth.

        First run self.calculate_impactForce() -> to find the force
        We don't run this everytime because it's time-consuming!

        :param azimuth: -> angle at which we're making 'virtual cut' to look at delaminations
        :param rmax: -> maximum r at which we do the delamination -> the larger this is the longer the analysis takes.
        :returns: a list of delamination lengths from bottom to top of laminate
        """
        rinterval = 0.01 # mm

        laminate = self.panel

        # Delamination lengths should be a list with length of the laminas list +1, each lamina has 2 sides!
        # However, the top and bottom should always have length zero as there is no delamination possible on the outside
        delaminationlengths = [0 for _ in range(len(laminate.laminas) + 1)]

        # we also define a 'last theta' -> angle of fibers at last lamina, (placeholder value = 0)
        lasttheta = 0

        for laminaindex, lamina in enumerate(laminate.laminas):
            # We need the interlaminar shear strength of the plies in the azimuth direction:
            Taucrit = lamina.Taucrit(azimuth)

            # check the shear stress at ply interfaces at different locations r:
            z0 = lamina.z0
            z1 = lamina.z1

            # Find the maximum shear stress:
            MaxTaurz0 = self.MaxTaurz(z0)
            MaxTaurz1 = self.MaxTaurz(z1)

            # Only do delamination analysis if it's even possible for a delamiination to occur:
            if MaxTaurz0[0] > Taucrit:
                bottomdelamination_length = self.calculate_delamination_length(rmax, z0, Taucrit, MaxTaurz0[1])
            else:
                bottomdelamination_length = 0

            # only add the delamination if the delamination is greater than previously found
            if bottomdelamination_length > delaminationlengths[laminaindex]:
                delaminationlengths[laminaindex] = bottomdelamination_length

            # However! if the current lamina has the same ply orientation as the previous (lamina below this), no delamination occurs!:
            # -> overwrite last delamination length
            if lamina.theta_ == lasttheta:
                delaminationlengths[laminaindex] = 0

            # Only do delamination analysis if it's even possible for a delamination to occur:
            if MaxTaurz1[0] > Taucrit:
                topdelamination_length = self.calculate_delamination_length(rmax, z1, Taucrit, MaxTaurz1[1])
            else:
                topdelamination_length = 0

            # only add the delamination if the delamination is greater than previously found
            if topdelamination_length > delaminationlengths[laminaindex + 1]:
                delaminationlengths[laminaindex + 1] = topdelamination_length

            lasttheta = lamina.theta_
        return delaminationlengths

    def calculate_delamination_length(self, rmax, z, Taucrit, r_maxTau):
        """
        Calculates delamination length, given a shear strength and a z coordinate
        :return:
        """
        # Create the array of data, Tau, r, integral of Tau as rows:
        stepsize = 0.005
        r_values = np.arange(0, rmax, stepsize)
        n = len(r_values)

        # array has 3 columns:
        # Column 1 is r
        # Column 2 is the function value
        # Column 3 is the integral from zero up to that point

        integral_value = 0
        integral_values = []
        for r in r_values:
            # Save Taurz
            Taur = self.Taurz(r, z)

            # Update integral and save
            # TODO: figure out wether to use cylindrical coordinates for this integral! it's a force equilibrium
            # the unit of the integral 'force' is in N per unit angle, as usually it would be in N per unit width if it
            # were a linear problem..
            integral_value += Taur * r * stepsize
            integral_values.append(integral_value)


        integral_interp = interp.interp1d(r_values, integral_values, kind='linear')

        # Function to find roots for
        def func_to_solve(r):
            return self.Taurz(r, z) - Taucrit

        # Find r values where Taurz equals Taucrit
        roots = []
        try:
            # Use root_scalar to find the roots
            result = opt.root_scalar(func_to_solve, bracket=[0, r_maxTau], method='brentq')
            if result.converged:
                roots.append(result.root)
        except ValueError:
            pass

        # It's known the function Taurz first goes above Taucrit and then dips back down again
        # We try to find the second root by searching beyond the first root
        if roots:
            try:
                result = opt.root_scalar(func_to_solve, bracket=[r_maxTau, rmax], method='brentq')
                if result.converged:
                    roots.append(result.root)
            except ValueError:
                pass

        # Get the integral values at the roots
        integral_at_roots = [integral_interp(root) for root in roots]

        # Now calculate the area (quasi volume) above the Taucrit and left of Taucrit:
        Areserve = (roots[0]**2) * Taucrit - integral_at_roots[0]
        Aexcess = (integral_at_roots[1] - integral_at_roots[0]) - (roots[1]**2-roots[0]**2) * Taucrit

        # If Aexcess is smaller than Areserve, then the delamination length is simply the last root:
        if Aexcess < Areserve:
            return roots[1]

        # If this is not the case, we have to calculate at which point the force equilibrium has been reached,
        # Beyond the last root!

        # -> find r at which the integral equals the area Taucrit * r:
        def AreaEq(r):
            return integral_interp(r) - Taucrit * r**2
        result = opt.root_scalar(AreaEq, bracket=[roots[1], rmax-stepsize], method='brentq')
        delaminationlength = result.root

        return delaminationlength

    def MaxTaurz(self, z):
        """
        returns the maximum value for Taurz given an Fimpact and z coordinate. This to determine wether
        any delamination happens at all
        :return:
        """

        # Interval in which you want to find the maximum for the first parameter
        bounds = (0, self.calculate_Rc())

        # Creating a lambda function to fix y and vary x
        fixed_y_func = lambda r: self.Taurz(r, z)

        # Minimizing the negative of the function to find the maximum
        result = opt.minimize_scalar(lambda x: -fixed_y_func(x), bounds=bounds, method='bounded')
        max_r = result.x

        MaxTaurz = self.Taurz(max_r, z)
        return MaxTaurz, max_r

    # Now we have to find direction in which damaged area is largest:
    def major_minor_axes(self):
        # ----------------------------------------------------------------------------
        # finding max delamination length:
        # ----------------------------------------------------------------------------
        angles = np.arange(0, 365, 5)
        lengths_angles = []
        for angle in tqdm(angles, desc = 'Delamination at all angles:'):
            delamination_lengths = self.delamination_analysis(angle, 15)
            maxdelamination_length = max(delamination_lengths)
            lengths_angles.append(maxdelamination_length)

        # Find the maximum value in lengths_angles
        max_length = max(lengths_angles)
        self.RadiusDelamination = max_length

        # Now find the angle at which the delamination is largest:
        max_index = lengths_angles.index(max_length)

        # Get the corresponding angle
        major_axis = angles[max_index]
        minor_axis = major_axis + 90

        # # Convert angles to radians for the polar plot
        # angles_radians = np.radians(angles)
        #
        # # Create the polar plot
        # plt.figure()
        # ax = plt.subplot(111, polar=True)
        # # ax.plot(angles_radians, lengths_angles)
        #
        # data = [self.delamination_analysis(angle, 15) for angle in angles]
        # # Plot each index of the returned lists
        # for i in range(len(data[0])):
        #     values = [d[i] for d in data]
        #     ax.plot(angles_radians, values, label=f'Index {i}')
        #
        # # Optionally, you can add labels and a title
        # ax.set_title('Delamination Lengths vs. Angles')
        # ax.set_xlabel('Angle (degrees)')
        # ax.set_ylabel('Delamination Length')
        #
        # # Show the plot
        # plt.show()

        # Return the major and minor axis directions
        return major_axis, minor_axis

    def generate_damaged_region(self,rmax):
        # We need to make zones, and one zone has a number of delaminations:
        delaminations = self.delamination_analysis(0,  rmax)

        # Now take out the zeros:
        delaminationlengths = [value for value in delaminations if value != 0.0]
        # Sort the list
        delaminationlengths = sorted(delaminationlengths)

        # Make the actual zones:
        zones = []
        for length in delaminationlengths:
            # In the delaminations list, zet equal to zero those lengths that are shorter than
            # current length:
            zonedelaminations = [value if value >= length else 0.0 for value in delaminations]
            # Now given these lengths generate the sublaminates:
            zone = self.generate_zone(zonedelaminations, length)
            zones.append(zone)

        self.damagedregion = DamagedRegion(zones)
        return self.damagedregion

    def calculate_CAI(self, knockdown, rmax):
        self.generate_damaged_region(rmax)
        E11 = self.panel.Ex
        E22 = self.panel.Ey
        v12 = self.panel.vxy
        G12 = self.panel.Gxy

        # Calculate the ratio between E moduli of laminate vs damaged region
        Eave = self.damagedregion.E_reduced * knockdown
        Elaminate = self.panel.Ex
        l = Elaminate / Eave

        # Find the size of the damage, we'll use the major axis in this case
        R = self.RadiusDelamination
        w = self.b

        num = 1 + (l + (1 - l * v12 ** 2 * (E22 / E11)) * np.sqrt(2 * (np.sqrt(E11 / E22) - v12) + (E11 / G12))) + (
                    (E11 / G12) - v12) * np.sqrt(E22 / E11)
        den = 1 + l * (l + (1 + np.sqrt(E22 / E11) * np.sqrt(2 * (np.sqrt(E11 / E22) - v12) + (E11 / G12)))) + (
                    (E11 / G12) - 2 * l * v12) * np.sqrt(E22 / E11) - (1 - l) ** 2 * v12 ** 2 * (E22 / E11)
        SCF1 = 1 - (1 - l) * num / den

        # Obtain the final stress concentration factor by compensating for the finite width:
        SCF2 = (2 + (1 - (2 * R / w) ** 3)) / (3 * (1 - (2 * R / w))) * SCF1
        return SCF2

    def generate_zone(self, delaminationlengths, length1):
        # make list of angles of the whole laminate:
        angles = []
        for lamina in self.panel.laminas:
            angles.append(lamina.theta_)

        # Remove the first element
        delaminationlengths = delaminationlengths[1:-1]
        delaminationlengths.append(1)

        # Now we have a list of angles which we just have to split at the delaminations:
        sublaminates_angleslists = []
        sublaminateangles = []
        for index, length in enumerate(delaminationlengths):
            if length == 0:
                sublaminateangles.append(angles[index])
            if length != 0:
                sublaminateangles.append(angles[index])
                sublaminates_angleslists.append(sublaminateangles)
                sublaminateangles = []

        sublaminates = []
        for sublaminateangles in sublaminates_angleslists:
            sublaminate = laminate_builder(sublaminateangles, False, False, 1)
            sublaminates.append(sublaminate)
        zone = Zone(sublaminates, length1, self.panel.h)
        return zone

    def plot_delamination(self, delaminationlengths, azimuth):
        Rc = self.calculate_Rc()
        fig, ax = plt.subplots()
        laminate = self.panel
        # Set x-limit beyond Rc
        xlim_value = max(delaminationlengths) * 2
        ax.set_xlim(0, xlim_value)

        # Plotting black horizontal lines for ply interfaces
        for i, lamina in enumerate(laminate.laminas):
            z0 = lamina.z0
            z1 = lamina.z1
            ax.plot([0, xlim_value], [z0, z0], color='black')
            ax.plot([0, xlim_value], [z1, z1], color='black')

        # Plotting blue lines for delamination lengths with increased thickness
        for i, length in enumerate(delaminationlengths):
            if i == 0 or i == len(delaminationlengths) - 1:
                continue  # Skip the top and bottom outer surfaces
            z = laminate.laminas[i - 1].z1  # Bottom of the current lamina
            ax.plot([0, length], [z, z], color='black', linewidth=4)

        # Plotting the vertical black dashed line at Rc
        ax.axvline(x=Rc, color='grey', linestyle='--', label='Rc')

        ax.set_xlabel('r (mm)')
        ax.set_ylabel('Ply Interface Position (z), (mm)')
        ax.set_title('Delamination Analysis at angle {}'.format(azimuth))
        plt.legend()
        plt.show()

    def plot_Taurz(self, z, rmax):
        r_values = np.linspace(0, rmax, 1000)
        Taurz_values = [self.Taurz(r, z) for r in r_values]
        plt.figure()
        plt.plot(r_values, Taurz_values, label=f'Taurz at z={z}')
        plt.axvline(x=self.calculate_Rc(), color='black', linestyle='--', label='Rc')
        plt.xlabel('r (mm)')
        plt.ylabel('Taurz')
        plt.title(f'Taurz as a function of r at z={z}')
        plt.legend()
        plt.show()

    def plot_force_equilibrium(self, z, rmax):
        r_values = np.linspace(0, rmax, 1000)
        Taurz_values = [self.Taurz(r, z) * r for r in r_values]
        Srz_values = [92 * r for r in r_values]
        plt.figure()
        plt.plot(r_values, Taurz_values, label=f'Taurz*r at z={z}')
        plt.plot(r_values, Srz_values, label=f'Srz*R at z={z}')
        plt.axvline(x=self.calculate_Rc(), color='black', linestyle='--', label='Rc')
        plt.xlabel('r (mm)')
        plt.ylabel('Taurz*r, Srz*r')
        plt.title(f'Integrands as a function of r at z={z}')
        plt.legend()
        plt.show()