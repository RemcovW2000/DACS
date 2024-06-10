import numpy as np
class Core:
    def __init__(self, thickness, coreproperties):
        self.thickness = thickness
        self.coreproperties = coreproperties

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
