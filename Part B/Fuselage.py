

class Fuselage:
    def __init__(self, diameter, n_joints, n_stringers):
        self.diameter        = diameter
        self.n_joints        = n_joints
        self.n_stringers     = n_stringers
        self.theta_joints    = self.spacing_joints()
        self.theta_stringers = self.spacing_stringers()
        # physical structural elements
        self.stringers       = []
        self.skins           = []
        # idealized structural elements
        self.booms           = []
        self.panels          = []
        # idealized 2nd moment of area
        self.Ixx             = 0

    def spacing_joints(self):
        if self.n_joints != 0:
            theta_joints  = 360/(2*self.n_joints)
            return theta_joints
        else:
            return None
        

