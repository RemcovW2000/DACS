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
Skin_2     = Skin.Skin_1
Skin_3     = Skin.Skin_1
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

failure_stringers, load_stringers_v, load_stringers_h, fpf_stringers_v            = [],[],[],[]
fpf_stringers_h, buckling_stringers, crippling_stringers_v, crippling_stringers_h = [],[],[],[]
failure_panels, load_panels_Nx, load_panels_Ns, fpf_panels, buckling_panels       = [],[],[],[],[]

data_stringers = {'Failure?':failure_stringers, 'Nx_v [N/mm]': load_stringers_v, 'Nx_h [N/mm]': load_stringers_h,  
                  'FPF_v': fpf_stringers_v, 'FPF_h': fpf_stringers_h,'Buckling':buckling_stringers,
                  'Crippling_v': crippling_stringers_v, 'Crippling_h': crippling_stringers_h}
data_panels    = {'Failure?':failure_panels, 'Nx [N/mm]': load_panels_Nx, 'Ns [N/mm]': load_panels_Ns, 
                  'FPF': fpf_panels, 'Buckling':buckling_panels}

for stringer in fuselage.stringers:
    stringer.FailureAnalysis()
    failure_stringers.append(stringer.Failure)
    load_stringers_v.append(np.round(stringer.Nxv,3))
    load_stringers_h.append(np.round(stringer.Nxh,3))
    fpf_stringers_v.append(np.round(stringer.FPFFIv,3))
    fpf_stringers_h.append(np.round(stringer.FPFFIh,3))
    buckling_stringers.append(np.round(stringer.BucklingFI,3))
    crippling_stringers_v.append(np.round(stringer.CripplingFIv,3))
    crippling_stringers_h.append(np.round(stringer.CripplingFIh,3))

df_stringers = pd.DataFrame(data_stringers)

for panel in fuselage.panels:
    panel.FailureAnalysis()
    failure_panels.append(panel.Failure)
    load_panels_Nx.append(np.round(panel.Nx,3))
    load_panels_Ns.append(np.round(panel.Ns,3))
    fpf_panels.append(np.round(panel.FPFFI,3))
    buckling_panels.append(np.round(panel.BucklingFI,3))

df_panels = pd.DataFrame(data_panels)
print('SUMMARY FAILURE ANALYSIS --------------------------------------------------------------------------------------------------------')
print('Stringers:')
print(df_stringers)
print('---------------------------------------------------------------------------------------------------------------------------------')
print('Panels:')
print(df_panels)
print('END------------------------------------------------------------------------------------------------------------------------------')


# ---------------------------------------------------------------------