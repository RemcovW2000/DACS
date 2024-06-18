import numpy as np
from Toolbox.Lamina import Lamina
from Toolbox import MP
from Toolbox.Laminate import Laminate, LaminateBuilder
import matplotlib.pyplot as plt

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
        Eindentation =0.4*k*(deltamax)**0.4
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
        # First run the following:
        # 1. self.impactforce()
        # We don't run this everytime because it's time consuming!
        rinterval = 0.001 # mm

        laminate = self.panel

        # Delamination lengths should be a list with length of the laminas list +1, each lamina has 2 sides!
        # However, the top and bottom should always have length zero as there is no delamination possible on the outside
        delaminationlengths = [0 for _ in range(len(laminate.laminas) + 1)]

        for laminaindex, lamina in enumerate(laminate.laminas):
            r = 0
            prev_topDL = True
            prev_botDL = True

            # We need the interlaminar shear strength of the plies in the azimuth direction:
            Taucrit = lamina.Taucrit(azimuth)

            while r < rmax:
                # check the shear stress at ply interfaces at different locations r:
                z0 = lamina.z0
                z1 = lamina.z1
                Tauz0 = self.Taurz(r, z0)
                Tauz1 = self.Taurz(r, z1)
                bottomdelamination, topdelamination = lamina.delaminationanalysis(Tauz0, Tauz1, azimuth)

                # When the function returns false for the firs time, the delamination has stopped at that r:
                if bottomdelamination != prev_botDL:
                    # print(Tauz0, 'at radius:', r)
                    bottomdelamination_length = r

                    # Now we'll add this to the delamination lengths list, if it is greater than the last delamination length
                    if bottomdelamination_length > delaminationlengths[laminaindex]:
                        delaminationlengths[laminaindex] = bottomdelamination_length
                    else:
                        pass

                if topdelamination != prev_topDL:
                    topdelamination_length = r

                    # Now we'll add this to the delamination lengths list, if it is greater than the last delamination length
                    if topdelamination_length > delaminationlengths[laminaindex + 1]:
                        delaminationlengths[laminaindex + 1] = topdelamination_length
                    else:
                        pass

                if not topdelamination and not bottomdelamination:
                    break

                prev_botDL = bottomdelamination
                prev_topDL = topdelamination

                r = r + rinterval
        return delaminationlengths

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
        r_values = np.linspace(2, rmax, 1000)
        Taurz_values = [self.Taurz(r, z) for r in r_values]

        plt.figure()
        plt.plot(r_values, Taurz_values, label=f'Taurz at z={z}')
        plt.axvline(x=self.calculate_Rc(), color='black', linestyle='--', label='Rc')
        plt.xlabel('r (mm)')
        plt.ylabel('Taurz')
        plt.title(f'Taurz as a function of r at z={z}')
        plt.legend()
        plt.show()


Laminate = LaminateBuilder([0, 90, 45, -45], True, True, 5)

Member = Member(Laminate)
impactforce = Member.impactforce(1000*1000, 150, 100)
deflection = Member.compute_deflection(impactforce, 150, 100, 150, 100)
print('Impact force: ', impactforce, 'N')
print('Deflection at this load: ', deflection, 'mm')

delamination_lengths = Member.DelaminationAnalysis(45, 20)
print(np.round(delamination_lengths, 2))

Member.plot_delamination(delamination_lengths)
Member.plot_Taurz(z=0, rmax=20)










