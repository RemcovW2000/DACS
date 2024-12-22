import numpy as np
import Stringer
import Skin

class Fuselage:
    def __init__(self, diameter, frame_spacing, n_joints, n_stringers, mass_frame = 50):
        self.diameter                  = diameter
        self.n_joints                  = n_joints
        self.n_stringers               = n_stringers
        self.theta_joints_array        = self.spacing_joints()
        self.theta_stringers_array     = self.spacing_stringers()
        self.theta_booms_array         = np.array([])
        # physical structural elements
        self.stringers       = []
        self.skins           = []
        self.frame_spacing   = frame_spacing
        # idealized structural elements
        self.booms           = []
        self.panels          = []
        # idealized 2nd moment of area
        self.Ixx             = 0
        self.mass_frame      = mass_frame # in grams
        self.AssignDiameter()

    def AssignDiameter(self):
        for i in self.skins:
            i.diameter = self.diameter
        return

    def spacing_joints(self): #NOTE: constant spacing for now
        if self.n_joints != 0:
            theta_joints       = 360/(2*self.n_joints)
            theta_joints_array = np.arange(theta_joints, 360, 2*theta_joints)
            return theta_joints_array
        else:
            return np.array([])
        
    def spacing_stringers(self): #NOTE: constant spacing for now
        if self.n_stringers != 0:
            theta_stringers       = 360/(2*self.n_stringers)
            theta_stringers_array = np.arange(theta_stringers, 360, 2*theta_stringers)
            custom_stringers_array = [18, 36, 144, 162, 198, 216, 324, 342]
            #theta_stringers_array = np.sort(np.concatenate((theta_stringers_array, custom_stringers_array)))
            #print(theta_stringers_array)
            return theta_stringers_array
        else:
            return np.array([])

    def Calculatempl(self):
        self.AssignDiameter()
        mpl = 0

        for i in self.stringers:
            mpl += i.Calculatempl()
        for i in self.skins:
            mpl += i.Calculatempl()

        print('frame spacing:', self.frame_spacing)
        print(mpl)
        self.mpl = (self.mass_frame + self.frame_spacing * mpl)/self.frame_spacing

        return self.mpl
    
    def CalculateNeutralAxis(self):
        tot_area  = 0
        static_moment = 0
        for stringer in self.stringers:
            tot_area             += stringer.area 
            static_moment        += stringer.area*np.sin(np.radians(stringer.location))*self.diameter/2 

        for skin in self.skins:
            tot_area      += skin.A
            static_moment += skin.static_moment

        self.neutral_axis = static_moment/tot_area

        return self.neutral_axis 