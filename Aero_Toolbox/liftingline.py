import numpy as np

class liftingline:
    def __init__(self, CLA, s, Cl_a):
        self.s = s      # semispan/halfspan of wing
        self.CLA = CLA  # average Cl times area of wing -> multiplied by dynamic pressure gives lift
        self.Cl_a = Cl_a # cl alpha slope

    def build_matrix(self, ntheta, nn):
        thetas = np.linspace(0.005*np.pi, 0.995*np.pi, ntheta)
        ns = [2 + 2*n for n in range(nn)]
        matrix = np.zeros((ntheta, nn))
        for i, theta in enumerate(thetas):
            for p, n in enumerate(ns):
                matrix[i, p] = self.left_side_eq(n, theta, self.chord_at(theta))
        self.matrix = matrix
        return matrix

    def left_side_eq(self, n, theta, c):
        s = self.s

        right_factor = self.Cl_a*np.sin(theta)/(8*s)
        left_val = np.sin(n*theta)*(np.sin(theta) + n*self.Cl_a*c/(8*s))
        value = left_val/right_factor
        return value

    def chord_at(self, theta):
        chord = 1
        return chord

r = Cl_a * alpha * c
r = Cl_a * (alphageom - alfai) * c
alphageom - alphai = r/(Cl_a * c)
alphai = alphageom - r/(Cl_a * c)
liftline = liftingline(10, 20, 2*np.pi)
print(liftline.build_matrix(10, 10))