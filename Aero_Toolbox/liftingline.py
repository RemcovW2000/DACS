import numpy as np

class liftingline:
    def __init__(self, ClA, s, Cl_a):
        self.s = s      # semispan/halfspan of wing
        self.CLA = CLA  # average Cl times area of wing -> multiplied by dynamic pressure gives lift
        self.Cl_a = Cl_a # cl alpha slope

    def build_matrix(self, ntheta, nn):
        thetas = np.linspace(0, 0.995*np.pi, ntheta)
        ns = [2*n for n in range(len(nn))]
        print(ns, thetas)
        return

    def left_side_eq(self, n, theta, c):
        right_factor = self.Cl_a*sin(theta)/(8*s)
        left_val = np.sin(n*theta)*(sin(theta) + n*self.Cl_a*c/(8*s))
        value = left_val/right_factor
        return value


liftline = liftingline(10, 20, 2*np.pi)
liftline.build_matrix(10, 10)