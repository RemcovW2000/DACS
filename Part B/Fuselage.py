import numpy as np

class Fuselage:
    def __init__(self, diameter, frame_spacing, n_joints, n_stringers):
        self.diameter        = diameter
        self.n_joints        = n_joints
        self.n_stringers     = n_stringers
        self.theta_joints    = self.spacing_joints()
        self.theta_stringers = self.spacing_stringers()
        # physical structural elements
        self.stringers       = []
        self.skins           = []
        self.frame_spacing   = frame_spacing
        # idealized structural elements
        self.booms           = []
        self.panels          = []
        # idealized 2nd moment of area
        self.Ixx             = 0

    def spacing_joints(self): #NOTE: constant spacing for now
        if self.n_joints != 0:
            theta_joints  = 360/(2*self.n_joints)
            theta_joints_array = np.arange(theta_joints, 360, 2*theta_joints)
            return theta_joints_array
        else:
            return []
        
    def spacing_stringers(self): #NOTE: constant spacing for now
        if self.n_stringers != 0:
            theta_stringers       = 360/(2*self.n_stringers)
            theta_stringers_array = np.arange(theta_stringers, 360, 2*theta_stringers)
            return theta_stringers_array
        else:
            return []

