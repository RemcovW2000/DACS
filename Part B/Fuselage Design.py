import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Fuselage import Fuselage
from Boom import Boom 
#from Skin import Skin 
from Stringer import Stringer
from Panel import Panel
from Structural_Idealization import Structural_Idealization
# ---------------------------------------------------------------------
# Remco:
import Stringer
import Skin
# Code
Stringer_1 = Stringer.TStringer_1
Stringer_2 = Stringer.TStringer_2
Stringer_3 = Stringer.TStringer_3
Skin_1     = Skin.Skin_1
Skin_2     = Skin.Skin_2
Skin_3     = Skin.Skin_3
# ---------------------------------------------------------------------

stringers = [Stringer_2, Stringer_3, Stringer_3, Stringer_3, Stringer_2, Stringer_1, Stringer_1, Stringer_1]
skins     = [Skin_1, Skin_2, Skin_3, Skin_1]
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
    buckling_stringers.append(stringer.BucklingFI)
    crippling_stringers.append((stringer.CripplingFIv, stringer.CripplingFIh))

df_stringers = pd.DataFrame(data_stringers)

for panel in fuselage.panels:
    panel.FailureAnalysis()
    failure_panels.append(panel.Failure)
    load_panels.append((panel.Nx, panel.Ns))
    fpf_panels.append(panel.FPFFI)
    buckling_panels.append(panel.BucklingFI)

df_panels = pd.DataFrame(data_panels)

print('Stringers:')
print(df_stringers)
print('-------------------------------------------------------------------------------')
print('Panels:')
print(df_panels)
print('-------------------------------------------------------------------------------')


# ---------------------------------------------------------------------