import numpy as np
from Toolbox.Lamina import Lamina
from Toolbox import MP
from Toolbox.Laminate import Laminate, LaminateBuilder
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp

class Member:
    def __init__(self, panel, loads = [0, 0, 0, 0, 0, 0], a = 300, b = 200):
        self.panel = panel # could be either laminate or sandwich, both should work
                            # Sandwich impact would not work
        self.loads = loads # total loads, not intensity!!!
        self.h = self.panel.h

        self.a = a # panel width
        self.b = b # panel depth

        self.R_impactor = 6.35
        self.E_impactor = 200000  # MPa
        self.v_impactor = 0.28

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
        max_m = 20
        max_n = 20

        result = 0
        for m in range(1, max_m + 1):
            for n in range(1, max_n + 1):
                result += self.deflection_single_term(F, m, n, xo, yo, x, y)
        return result

# ----------------------------------------------------------------------------------------------------------------------
# DACS II impact damage:
    # Find the force:
# ----------------------------------------------------------------------------------------------------------------------
    def Edeflection(self, F, xo, yo):
        # For the deflection energy we analyse the work done at location of impact so
        # we need the deflection and force at place of impact:
        x = xo
        y = yo

        wmax = self.compute_deflection(F, xo, yo, x, y)
        Edeflection = F*wmax/2
        return Edeflection

    def Eindentation(self, F):
        k = self.k_stiffness(self.R_impactor)
        deltamax = self.calculate_deltamax(F, k)
        Eindentation =0.4*k*(deltamax)**2.5
        return Eindentation

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

        ABD = self.panel.ABD_matrix
        A11 = ABD[0,0]
        A22 = ABD[1,1]
        A12 = ABD[0,1]

        # Ask how to calculate this:
        Gzr = 4500 # PLACEHOLDER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        term1 = np.sqrt(A22)
        inner_term1 = np.sqrt(A11 * A22) + Gzr
        inner_term2 = A12 + Gzr
        numerator = term1 * np.sqrt(inner_term1 ** 2 - inner_term2 ** 2)

        # Calculate components of the denominator
        denominator = 2 * np.pi * np.sqrt(Gzr) * (A11 * A22 - A12 ** 2)
        K2 = numerator / denominator

        K2 = (1-0.3**2)/11.2e3
        k = 4*np.sqrt(R)/(3*np.pi*(K1 + K2))
        return k

    def Eimpact_estimate(self, F, xo, yo):
        return self.Eindentation(F) + self.Edeflection(F, xo, yo)

    def impactforce(self, Eimpact, xo, yo, tol=1e-4, max_iter=1000):
        # This function sets the attribute Fimpact and also returns the force

        # Let's also save the impact energy for which
        self.Eimpact = Eimpact

        F = 1000 # Initial guess for F
        for i in range(max_iter):
            E_est = self.Eimpact_estimate(F, xo, yo)
            diff = E_est - Eimpact

            if abs(diff) < tol:
                # Set impact force attribute for later use:
                self.Fimpact = F

                # Set delta max attribute for later use:
                self.deltamax = self.calculate_deltamax(F, self.k_stiffness(self.R_impactor))
                return F  # Found the solution within tolerance

            # Newton-Raphson update step
            dE_dF = (self.Eimpact_estimate(F + tol, xo, yo) - E_est) / tol  # Numerical derivative
            F = F - diff / dE_dF

            if F < 0:
                F = tol  # Ensure F stays positive, you can choose a small positive number

        raise ValueError("Failed to converge to a solution")

    def Taurz(self, r, z):
        # First run the following:
        # 2. self.impactforce()
        # We don't run this everytime because it's time consuming!
        Ftotal = self.Fimpact
        Rc = self.calculate_Rc()

        if r < Rc:
            if r == 0:
                r = 1e-10
            a = 1/Rc**2
            term1 = 3*a*Ftotal
            term2 = 1/(3*a)
            numerator = (1-a*r**2)**1.5
            denominator = 3*a
            Fr = term1*(term2 -numerator/denominator)
            Tauavg = Fr/(2*np.pi*r*self.h)
        elif r >= Rc:
            Tauavg = Ftotal/(2*np.pi*r*self.h)

        Taumax = 1.5*Tauavg

        # Now use Tauavg to find the
        Taurz = Taumax*(1-((2*z)**2)/self.h**2)
        return Taurz

    def DelaminationAnalysis(self, azimuth, rmax):
        """
        Checks delaminations at a specific angle azimuth.

        :param azimuth: -> angle at which we're making 'virtual cut' to look at delaminations
        :param rmax: -> maximum r at which we do the delamination -> the larger this is the longer the analysis takes.
        :returns: a list of delamination lengths from bottom to top of laminate
        """
        # First run the following:
        # 1. self.impactforce()
        # We don't run this everytime because it's time-consuming!
        rinterval = 0.001 # mm

        laminate = self.panel

        # Delamination lengths should be a list with length of the laminas list +1, each lamina has 2 sides!
        # However, the top and bottom should always have length zero as there is no delamination possible on the outside
        delaminationlengths = [0 for _ in range(len(laminate.laminas) + 1)]

        # we also define a 'last theta' -> angle of fibers at last lamina, (placeholder value = 0)
        lasttheta = 0

        for laminaindex, lamina in enumerate(laminate.laminas):
            r = 0
            prev_topDL = True
            prev_botDL = True

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
                bottomdelamination_length = Member.calculate_delamination_length(rmax, z0, Taucrit, MaxTaurz0[1])
            else:
                bottomdelamination_length = 0

            # only add the delamination if the delamination is greater than previously found
            if bottomdelamination_length > delaminationlengths[laminaindex]:
                delaminationlengths[laminaindex] = bottomdelamination_length

            # However! if the current lamina has the same ply orientation as the previous (lamina below this), no delamination occurs!:
            # -> overwrite last delamination length
            if lamina.theta == lasttheta:
                delaminationlengths[laminaindex] = 0

            # Only do delamination analysis if it's even possible for a delamination to occur:
            if MaxTaurz1[0] > Taucrit:
                topdelamination_length = Member.calculate_delamination_length(rmax, z1, Taucrit, MaxTaurz1[1])
            else:
                topdelamination_length = 0

            # only add the delamination if the delamination is greater than previously found
            if topdelamination_length > delaminationlengths[laminaindex + 1]:
                delaminationlengths[laminaindex + 1] = topdelamination_length
        return delaminationlengths

    def calculate_delamination_length(self, rmax, z, Taucrit, r_maxTau):
        """
        Calculates delamination length, given a shear strength and a z coordinate
        :return:
        """
        # Create the array of data, Tau, r, integral of Tau as rows:
        stepsize = 0.01
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
            integral_value += Taur * stepsize
            integral_values.append(integral_value)


        integral_interp = interp.interp1d(r_values, integral_values, kind='linear')

        # Function to find roots for
        def func_to_solve(r):
            return self.Taurz(r, z) - Taucrit

        # Find r values where Taurz equals specific_value
        roots = []
        try:
            # Use root_scalar to find the roots
            result = opt.root_scalar(func_to_solve, bracket=[0, r_maxTau], method='brentq')
            if result.converged:
                roots.append(result.root)
        except ValueError:
            pass

        # It's known the function Taurz first goes above the specific value and then dips back down again
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

        # Now calculate the area above the Taucrit and left of Taucrit:
        Areserve = roots[0]*Taucrit - integral_at_roots[0]
        Aexcess = (integral_at_roots[1] - integral_at_roots[0]) - (roots[1]-roots[0]) * Taucrit

        # If Aexcess is smaller than Areserve, then the delamination length is simply the last root:
        if Aexcess < Areserve:
            return roots[1]

        # If this is not the case, we have to calculate at which point the force equilibrium has been reached,
        # Beyond the last root!

        # -> find r at which the integral equals the area Taucrit * r:
        def AreaEq(r):
            return integral_interp(r) - Taucrit * r
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
    def Major_Minor_axes(self):
        angles = np.arange(0, 360, 5)
        lengths_angles = []
        for angle in angles:
            delamination_lengths = self.DelaminationAnalysis(angle, 50)
            maxdelamination_length = max(delamination_lengths)
            print(maxdelamination_length)
            lengths_angles.append(maxdelamination_length)

        # Find the maximum value in lengths_angles
        max_length = max(lengths_angles)

        # Now find the angle at which the delamination is largest:
        max_index = lengths_angles.index(max_length)

        # Get the corresponding angle
        major_axis = angles[max_index]
        minor_axis = major_axis + 90

        # Convert angles to radians for the polar plot
        angles_radians = np.radians(angles)

        # Create the polar plot
        plt.figure()
        ax = plt.subplot(111, polar=True)
        ax.plot(angles_radians, lengths_angles)

        # Optionally, you can add labels and a title
        ax.set_title('Delamination Lengths vs. Angles')
        ax.set_xlabel('Angle (degrees)')
        ax.set_ylabel('Delamination Length')

        # Show the plot
        plt.show()

        # Return the major and minor axis directions
        return major_axis, minor_axis

    def SublaminateBuilder(self):
        # First make a list of all the lamina angles:
        laminaangleslist = []
        for lamina in self.panel:
            laminaangleslist.append(lamina.theta)



        Sublaminates = []
        return Sublaminates

    def plot_delamination(self, delaminationlengths):
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
            ax.plot([0, length], [z, z], color='blue', linewidth=2)

        # Plotting the vertical black dashed line at Rc
        ax.axvline(x=Rc, color='black', linestyle='--', label='Rc')

        ax.set_xlabel('Delamination Length (mm)')
        ax.set_ylabel('Ply Interface Position (z)')
        ax.set_title('Delamination Analysis')
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


Laminate = LaminateBuilder([0, 90, 0, 0, 45, -45], True, True, 3)

Member = Member(Laminate)
impactforce = Member.impactforce(1000*30, 150, 100)
deflection = Member.compute_deflection(impactforce, 150, 100, 150, 100)
print('Impact force: ', impactforce, 'N')
print('Deflection at this load: ', deflection, 'mm')

delamination_lengths = Member.DelaminationAnalysis(45, 20)
print(np.round(delamination_lengths, 2))

Member.plot_delamination(delamination_lengths)
Member.plot_Taurz(z=0, rmax=20)

print(Member.Major_Minor_axes())










