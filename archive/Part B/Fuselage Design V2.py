import numpy as np
import copy
import matplotlib.pyplot as plt
import pandas as pd
#from Skin import Skin
from Structural_Idealization_V2 import Structural_Idealization_V2
# ---------------------------------------------------------------------
# Remco:
import Stringer
import Skin
# Code
Stringer_shear             = Stringer.TStringer_shear
Stringer_compression       = Stringer.TStringer_compression
Stringer_tension           = Stringer.TStringer_tension
Skin_compression     = Skin.Skin_compression
Skin_tension         = Skin.Skin_tension
Skin_shear           = Skin.Skin_shear
# ---------------------------------------------------------------------

stringers = [] # [copy.deepcopy(Stringer_1), copy.deepcopy(Stringer_1), copy.deepcopy(Stringer_1), copy.deepcopy(Stringer_1), copy.deepcopy(Stringer_1)]

n_custom_stringers = 0
n_stringers        = 0 + n_custom_stringers

while len(stringers) < n_stringers:

    if len(stringers) <= n_stringers/8:
        stringer = copy.deepcopy(Stringer_shear)
    elif n_stringers/8 < len(stringers) <= n_stringers*3/8:
        stringer = copy.deepcopy(Stringer_compression)
    elif n_stringers*3/8 < len(stringers) <= n_stringers*5/8:
        stringer = copy.deepcopy(Stringer_shear)
    elif n_stringers*5/8 < len(stringers) <= n_stringers*7/8:
        stringer = copy.deepcopy(Stringer_tension)
    elif n_stringers*7/8 < len(stringers) <= n_stringers:
        stringer = copy.deepcopy(Stringer_shear)

    stringers.append(stringer)

skins = [copy.deepcopy(Skin_compression), copy.deepcopy(Skin_shear), copy.deepcopy(Skin_tension), copy.deepcopy(Skin_shear)]

skinsmono = [copy.deepcopy(Skin_shear), copy.deepcopy(Skin_compression), copy.deepcopy(Skin_compression),
         copy.deepcopy(Skin_compression), copy.deepcopy(Skin_compression), copy.deepcopy(Skin_compression), copy.deepcopy(Skin_shear), copy.deepcopy(Skin_shear),
         copy.deepcopy(Skin_shear), copy.deepcopy(Skin_tension), copy.deepcopy(Skin_tension), copy.deepcopy(Skin_tension), copy.deepcopy(Skin_tension), copy.deepcopy(Skin_tension),
        copy.deepcopy(Skin_shear), copy.deepcopy(Skin_shear)]


# ---------------------------------------------------------------------
# fuselage parameters:
diameter      = 6e3 # [mm]
frame_spacing = 0.5e3 # [mm]
# load case
Vy = -1.5e6 # [N]
Mx = -15e9  # [Nmm]
mass_frame = 31.32 # kg

fuselage = Structural_Idealization_V2(Mx, Vy, diameter, frame_spacing, stringers, skinsmono, n_custom_stringers)
fuselage.mass_frame = mass_frame*1000 # weight of a frame in GRAMS!!!
fuselage.Calculatempl()


failure_stringers, load_stringers_v, load_stringers_h, fpf_stringers_v                                 = [],[],[],[]
fpf_stringers_h, buckling_stringers, crippling_stringers_v, crippling_stringers_h                      = [],[],[],[]
failure_panels, load_panels_Nx, load_panels_Ns, fpf_panels, buckling_panels, monolithic_buckling       = [],[],[],[],[],[]
stringer_locations, panel_locations                                                                    = [], []


data_stringers = {'Failure?':failure_stringers, 'Nx_v [N/mm]': load_stringers_v, 'Nx_h [N/mm]': load_stringers_h,
                  'FPF_v': fpf_stringers_v, 'FPF_h': fpf_stringers_h,'Buckling':buckling_stringers,
                  'Crippling_v': crippling_stringers_v, 'Crippling_h': crippling_stringers_h, 'Location': stringer_locations}
data_panels    = {'Failure?':failure_panels, 'Nx [N/mm]': load_panels_Nx, 'Ns [N/mm]': load_panels_Ns, 
                  'FPF': fpf_panels, 'Bu':buckling_panels, 'Mb': monolithic_buckling, 'Starts @': panel_locations}

data_fuselage  = {'MPL': [fuselage.mpl]}

df_fuselage    = pd.DataFrame(data_fuselage)

for stringer in fuselage.stringers:
    stringer.failure_analysis()
    f = 0
    if stringer.Failure == True:
        f += 1
        failure_stringers.append(stringer.Failure)
        load_stringers_v.append(np.round(stringer.Nxv,3))
        load_stringers_h.append(np.round(stringer.Nxh,3))
        fpf_stringers_v.append(np.round(stringer.FPFFIv,3))
        fpf_stringers_h.append(np.round(stringer.FPFFIh,3))
        buckling_stringers.append(np.round(stringer.BucklingFI,3))
        crippling_stringers_v.append(np.round(stringer.CripplingFIv,3))
        crippling_stringers_h.append(np.round(stringer.CripplingFIh,3))
        stringer_locations.append(str(stringer.location) + '\u00B0')
        
    if f == 0:
        if stringer.Failure == False:
            failure_stringers.append(stringer.Failure)
            load_stringers_v.append(np.round(stringer.Nxv,3))
            load_stringers_h.append(np.round(stringer.Nxh,3))
            fpf_stringers_v.append(np.round(stringer.FPFFIv,3))
            fpf_stringers_h.append(np.round(stringer.FPFFIh,3))
            buckling_stringers.append(np.round(stringer.BucklingFI,3))
            crippling_stringers_v.append(np.round(stringer.CripplingFIv,3))
            crippling_stringers_h.append(np.round(stringer.CripplingFIh,3))
            stringer_locations.append(str(stringer.location) + '\u00B0')
        


df_stringers = pd.DataFrame(data_stringers)

for panel in fuselage.panels:
    panel.failure_analysis()
    f = 0
    if panel.Failure == True:
        f += 1
        failure_panels.append(panel.Failure)
        load_panels_Nx.append(np.round(panel.Nx,3))
        load_panels_Ns.append(np.round(panel.Ns,3))
        fpf_panels.append(np.round(panel.FPFFI,3))
        buckling_panels.append(np.round(panel.BucklingFI,3))
        monolithic_buckling.append(np.round(panel.MonolithicBucklingFI, 3))
        panel_locations.append(str(panel.start) + '\u00B0')
    if f == 0:
        if  panel.Failure == False:
                failure_panels.append(panel.Failure)
                load_panels_Nx.append(np.round(panel.Nx,3))
                load_panels_Ns.append(np.round(panel.Ns,3))
                fpf_panels.append(np.round(panel.FPFFI,3))
                buckling_panels.append(np.round(panel.BucklingFI,3))
                monolithic_buckling.append(np.round(panel.MonolithicBucklingFI, 3))
                panel_locations.append(str(panel.start) + '\u00B0')


df_panels = pd.DataFrame(data_panels)

print('SUMMARY FAILURE ANALYSIS --------------------------------------------------------------------------------------------------------')
print('Stringers:')
print(df_stringers)
print('---------------------------------------------------------------------------------------------------------------------------------')
print('Panels:')
print(df_panels)
print('---------------------------------------------------------------------------------------------------------------------------------')
print('Fuselage:')
print(df_fuselage)

def show_structural_elements(ax, theta, color, label, marker):
    for angle in theta:
        x_dot = (diameter/2) * np.cos(np.radians(angle))
        y_dot = (diameter/2) * np.sin(np.radians(angle))
        if angle == theta[0]:
            ax.plot(x_dot, y_dot, marker, markersize = 13, color = color, label = label)
        else:
            ax.plot(x_dot, y_dot, marker,  markersize = 13, color = color)
        offset = 100
        ax.text(x_dot + 3*offset, y_dot + offset, f'{angle:.2f}°', fontsize=8, verticalalignment='bottom', horizontalalignment = 'left')

def show_fuselage(theta_s_array, theta_j_array, theta_b_array, y_bar):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle('FUSELAGE IDEALIZATION')
    axes = [ax1, ax2]
    # Generate circle points
    theta_circle = np.linspace(0, 2*np.pi, 100)
    x_circle = diameter/2 * np.cos(theta_circle)
    y_circle = diameter/2 * np.sin(theta_circle)

    # Plot circle
    for ax in axes:
        if ax == ax1:
            label = 'skin'
        else:
            label = 'panel'
        ax.plot(x_circle, y_circle, 'k-', label=label)

    # Plot Physical Structure
    if theta_s_array.size > 0:
        show_structural_elements(ax1, theta_s_array, 'r', label = 'stringers', marker='o')
    show_structural_elements(ax1, theta_j_array, 'b', label = 'joints', marker='x')
    # Plot Idealized Structure
    show_structural_elements(ax2, theta_b_array , 'gray', label = 'booms', marker='o')

    for ax in axes:
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.axis('equal')
        ax.grid(True)
        ax.legend()
        ax.axhline(y = y_bar, color='r', linestyle='--')
        ax.invert_xaxis()
        ax.invert_yaxis()
    ax1.set_title('PHYSICAL FUSELAGE')
    ax2.set_title('IDEALIZED FUSELAGE')

    plt.show()

print('y_bar', fuselage.neutral_axis)
show_fuselage(fuselage.theta_stringers_array, fuselage.theta_joints_array, fuselage.theta_booms_array, fuselage.neutral_axis)


