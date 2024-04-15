import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Fuselage import Fuselage
from Boom import Boom 
from Skin import Skin 
from Stringer import Stringer
from Panel import Panel
from Structural_Idealization import Structural_Idealization
# ---------------------------------------------------------------------
# Remco:
import Stringer
# Code
Stringer_a = Stringer.TStringer_a
Stringer_b = Stringer.TStringer_b
Stringer_c = Stringer.TStringer_c

# ---------------------------------------------------------------------

stringers = []
skins     = []
# ---------------------------------------------------------------------
# Killian:
# Code
# fuselage parameters:
diameter      = 6e3 # [mm]
frame_spacing = 2e3 # [mm]
# load case
Vy = -1.5e6 # [N]
Mx = -15e9  # [Nmm]

fuselage = Structural_Idealization(Mx, Vy, diameter, frame_spacing, stringers, skins)

failure_stringers, load_stringers, fpf_stringers, buckling_stringers, crippling_stringers = [],[],[],[],[]
failure_panels, load_panels, fpf_panels, buckling_panels                                  = [],[],[],[] 

data_stringers = {'Failure?':failure_stringers, 'Nx (v,h) [N/mm]': load_stringers, 
                  'FPF (v,h)': fpf_stringers, 'Buckling (v,h)':buckling_stringers, 
                  'Crippling (v,h)': crippling_stringers}
data_panels    = {'Failure?':failure_panels, 'Nx [N/mm]': load_panels, 
                  'FPF': fpf_panels, 'Buckling':buckling_panels}

for stringer in fuselage.stringers:
    stringer.FailureAnalysis()
    failure_stringers.append(stringer.Failure)
    load_stringers.append((stringer.Nxv, stringer.Nxh))
    fpf_stringers.append((stringer.FPFFIv, stringer.FPFFIh)) 
    buckling_stringers.append((stringer.BucklingFIv, stringer.BucklingFIh))
    crippling_stringers.append((stringer.CripplingFIv, stringer.CripplingFIh))

for panel in fuselage.panels:
    panel.FailureAnalysis()
    failure_panels.apppend(panel.Failure)
    load_panels.append((panel.Nx, panel.Ns))
    fpf_panels.append(panel.FPFFI)
    buckling_panels.append(panel.BucklingFI)






# ---------------------------------------------------------------------