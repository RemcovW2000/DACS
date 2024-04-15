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

for stringer in fuselage.stringers:
    stringer.FailureAnalysis()
    failure_stringers.append(stringer.Failure)
    load_stringers.append((stringer.Nxh, stringer.Nxv))
    fpf_stringers.append(stringer.FPFFI) 
    buckling_stringers.append(stringer.BucklingFI)
    crippling_stringers.append(stringer.CripplingFI)

for panel in fuselage.panels:
    panel.FailureAnalysis()
    failure_panels.apppend(panel.Failure)
    load_stringers.append(stringer)


# ---------------------------------------------------------------------