import numpy as np
class Core:
    def __init__(self, h, coreproperties):
        self.h = h
        self.coreproperties = coreproperties

        if self.coreproperties['G']:
            self.G = self.coreproperties['G']
        else:
            self.G = None

    def Gxbarz(self, theta):
        Gxz = self.coreproperties['Gxz']
        Gyz = self.coreproperties['Gyz']
        sin_theta_squared = np.sin(theta) ** 2
        cos_theta_squared = np.cos(theta) ** 2
        Gxbarz = sin_theta_squared * Gyz + cos_theta_squared * Gxz
        return Gxbarz

    def Gybarz(self, theta):
        Gxz = self.coreproperties['Gxz']
        Gyz = self.coreproperties['Gyz']
        sin_theta_squared = np.sin(theta) ** 2
        cos_theta_squared = np.cos(theta) ** 2
        Gybarz = cos_theta_squared * Gyz + sin_theta_squared * Gxz
        return Gybarz

    def calculate_weight_per_A(self):
        '''
        Calculates the weight of the core material per unit area based on the density and thickness. Does not take into
        account extra weight per unit area of the core material due to soaking of resin into partially open cells.

        :return:
        '''
        return self.h * self.coreproperties['rho']
