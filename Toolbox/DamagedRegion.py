import numpy as np

# Define a class for a damage zone, definition as given in the project description

class Zone:
    def __init__(self, sublaminates, width, height):
        self.sublaminates = sublaminates # idea is to receive sublaminates in a list
                                         # which are Laminate objects, such that Ex can 
                                         # be retrieved
                                         # TODO: check whether Ex non-symmetric laminate
                                         # is implemented
        self.width        = width        # width equals delamination size that defines
                                         # this corresponding zone
        self.height       = height       # total thickness laminate
        self.E_parallel   = self.get_E_parallel()

        self.E_equivalent = 0            # initialized to be zero, value depends on properties
                                         # of enclosed zones

    def get_E_parallel(self): 
        # calculates the E modulus through the thickness of a zone
        # in which the sublaminates acts as springs in parallel (strain compatibility)
    
        E_parallel = 0
        for sublaminate in self.sublaminates:
            E_parallel += sublaminate.h * sublaminate.Ex # sum all sublaminate contributions

        E_parallel = E_parallel/self.height # divide by total thickness
        return E_parallel
    

class DamagedRegion:
    def __init__(self, zones):
        self.zones     = zones 
        self.E_reduced = self.get_reduced_stiffness()

    def get_reduced_stiffness(self):
        for i, zone in enumerate(self.zones):
            if i == 0:
                zone.E_equivalent = zone.E_parallel # there is no zone enclosed in zone 1
            else:
                E_eq_enclosed      = self.zones[i-1].E_equivalent # E equivalent of enclosed zones
                d_inner            = self.zones[i-1].width        # delamination length defining enclosed zone
                d_outer            = zone.width                   # delamination length defining current zone
                zone.E_equivalent  = E_eq_enclosed * zone.E_parallel / (E_eq_enclosed * (1 - d_inner/d_outer) + zone.E_parallel*(d_inner/d_outer))
        E_reduced = self.zones[-1].E_equivalent

        return E_reduced






    




        