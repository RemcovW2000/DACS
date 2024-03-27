from ABD_Killian import Material, Lamina, Layup, Laminate, GenerateLaminate, ABD_Matrix, T
from StressAnalysis_Killian import Stress
from FailureCriteria_Killian import FibreFailure, MatrixFailure, MaxStressFailure
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
# print settings
np.set_printoptions(linewidth=300, precision = 3)

# regulate what questions are active
q_1a = False
q_1b = False
q_2a = True

'example analysis below'
# specify lay-up(including symmetry and thicknesses), materials
# lay_up    = Layup([0,45,-45,90], 1,  [0.2]) 
# material  = Material('UD CFRP', 140e3, 10e3, 5e3, 0.3, None, [1.5e3, 1.2e3, 0.05e3, 0.25e3, 0.07e3], 0.2, 230e3) # MPa, [-]
# print('lay_up:', lay_up.stack)

# # initialize lamina and laminate
# lamina   = Lamina(material, None, None, None, None, None, None, None)
# laminate = Laminate('cross-ply', lay_up.stack, None, np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3)), np.empty((6,6)), np.empty((6,6)), 0, 0, 0, 0, 0, np.empty((3,len(lay_up.stack))))

# # assign lamina orientations and thicknesses based upon specified lay-up 
# GenerateLaminate(lamina, lay_up, laminate) #, laminae)
# # calculate ABD
# ABD_Matrix(lay_up, laminate)

# # initialize force vector F, containing force and moment intensities
# F = np.ones((6,1)) # (Nx, Ny, Ns, Mx, My, Ms)
# # apply loadcase
# F[0] = 10 # N/mm -> results in stresses in MPa (N/mm^2)

# Stress(laminate, F)

# #print(laminate.stresses)

'1. Calculation of the ABD matrix and stress analysis (25%)'
'1a'
if q_1a:
    n = 1 # total number of plies = 2*n*5
    t_1a = 0.125 # wasn't given in assignment, assume this value and continue
    # initialize material and a lamina
    #material_1a  = Material('UD material', 145.3e3, 8.5e3, 4.58e3, 0.31, None, [1932, 0, 108, 0, 132.8], 0.2, 230e3) # MPa, [-]
    #print(material_1a.v21)
    'TEST TO IMPLEMENT UNCERTAINTIES'
    E1 = ufloat(145.3e3, 3.28e3)
    E2 = ufloat(8.5e3, 1.28e3)
    v12 = ufloat(0.31, 0.018)
    material_1a  = Material('UD material', E1, E2, 4.58e3, v12, None, [1932, 0, 108, 0, 132.8], 0.2, 230e3) # MPa, [-]
    lamina_1a    = Lamina(material_1a, None, None, None, None, None, None, None)

    theta_range = np.arange(0,91,1)
    Ex, Ey, Gxy, vxy, vyx = [],[],[],[],[]
    for theta in theta_range:
        base_stack = [15, theta, -theta, 75, 75]
        half_stack = base_stack
        for i in np.arange(1,n):
            half_stack += base_stack 
        lay_up_1a = Layup(half_stack, 1, [t_1a])
        # initialize laminate
        laminate_1a = Laminate('1a', lay_up_1a.stack, None, np.zeros((3,3) ,dtype='O'), np.zeros((3,3),dtype='O'), np.zeros((3,3),dtype='O'), np.empty((6,6)), np.empty((6,6)), 0, 0, 0, 0, 0, np.empty((3,len(lay_up_1a.stack))) ,np.empty((3,len(lay_up_1a.stack))))
        # assign lamina orientations and thicknesses based upon specified lay-up 
        GenerateLaminate(lamina_1a, lay_up_1a, laminate_1a) #, laminae)
        # calculate ABD and equivalent properties
        ABD_Matrix(lay_up_1a, laminate_1a)
        Ex.append(laminate_1a.Ex) # ADAPT FOR CHANGE OF DATA STORAGE
        Ey.append(laminate_1a.Ey)
        Gxy.append(laminate_1a.Gxy)
        vxy.append(laminate_1a.vxy)
        vyx.append(laminate_1a.vyx)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    string = ']'
    'Upper bound'
    ax1.plot(theta_range, [(i.n + i.s) / 1000 for i in Ex], color='r')
    ax1.plot(theta_range, [(i.n + i.s) / 1000 for i in Ey], color='r')
    ax1.plot(theta_range, [(i.n + i.s) / 1000 for i in Gxy], color='r')
    'Mean'
    ax1.plot(theta_range, [i.n / 1000 for i in Ex], label = 'Ex')
    ax1.plot(theta_range, [i.n / 1000 for i in Ey], label = 'Ey')
    ax1.plot(theta_range, [i.n / 1000 for i in Gxy], label = 'Gxy')
    'Lower bound'
    ax1.plot(theta_range, [(i.n - i.s) / 1000 for i in Ex],  color='r')
    ax1.plot(theta_range, [(i.n - i.s) / 1000 for i in Ey],  color='r')
    ax1.plot(theta_range, [(i.n - i.s) / 1000 for i in Gxy], color='r')
    'Labelling etc'
    ax1.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_{}'.format(2*n))
    ax1.set_ylabel('GPa')
    ax1.grid()
    ax1.legend()

    'Upper bound'
    ax2.plot(theta_range, [i.n + i.s for i in vxy], color='r')
    ax2.plot(theta_range, [i.n + i.s for i in vyx], color='r')
    'Mean'
    ax2.plot(theta_range, [i.n  for i in vxy], label = 'vxy')
    ax2.plot(theta_range, [i.n  for i in vyx], label = 'vyx')
    'Lower bound'
    ax2.plot(theta_range, [i.n - i.s  for i in vxy], color='r')
    ax2.plot(theta_range, [i.n - i.s  for i in vyx], color='r')

    ax2.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_{}'.format(2*n))
    ax2.set_ylabel('[-]')
    ax2.grid()
    ax2.legend()
    plt.show()

'1b' '''STILL NEED TO CONVERT STRESSES TO PRINCIPAL COO SYSTEM'''
if q_1b:
    # initialize material and a lamina
    material_1b  = Material('UD material', 145.3e3, 8.5e3, 4.58e3, 0.31, None, [1932, 0, 108, 0, 132.8], 0.2, 230e3) # MPa, [-]
    #material_1b  = Material('UD material', 140, 10, 5, 0.3, None, [1932, 0, 108, 0, 132.8], 0.2, 230e3) # MPa, [-]
    lamina_1b    = Lamina(material_1b, None, None, None, None, None, None, None, 0)
    #
    lay_up_1b = Layup([0,0,90,30,90], 0, [0.125])
    # initialize laminate
    laminate_1b = Laminate('1b', lay_up_1b.stack, None, np.zeros((3,3) ,dtype='O'), np.zeros((3,3),dtype='O'), np.zeros((3,3),dtype='O'), np.empty((6,6),dtype='O'), np.empty((6,6),dtype='O'), 0, 0, 0, 0, 0, np.empty((3,len(lay_up_1b.stack))), np.empty((3,len(lay_up_1b.stack))))
    # assign lamina orientations and thicknesses based upon specified lay-up 
    GenerateLaminate(lamina_1b, lay_up_1b, laminate_1b)
    # calculate ABD and equivalent properties
    ABD_Matrix(lay_up_1b, laminate_1b)
    # initialize force vector F, containing force and moment intensities
    F_1b = np.zeros((6,1)) # (Nx, Ny, Ns, Mx, My, Ms)
   
    # apply loadcase
    F_1b[0] = 0.2e-1 # N/mm -> results in stresses in MPa (N/mm^2)
    #F_1b[0]  = 100 # N/mm
    F_1b[1] = 1.8e1
    F_1b[3] = 18 #e3
    # calculate strains and stresses
    Stress(laminate_1b, F_1b)
    # results
    print(laminate_1b.ABD)
    print(laminate_1b.stresses) # MPa
    print(laminate_1b.strains)  # microstrains 

    t_coo = []
    for lamina in laminate_1b.laminae:
        t_coo.append(max([lamina.z0, lamina.z1], key=abs))

    print(t_coo)

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3) #, figsize=(10, 5))

    axis = [ax1, ax2, ax3, ax4, ax5, ax6]

    ax1.plot(laminate_1b.strains[0][:], t_coo)
    ax1.set_ylabel('z (mm)')
    ax1.set_xlabel('e_11 [microstrain]')
    ax2.plot(laminate_1b.strains[1][:], t_coo)
    ax2.set_xlabel('e_22 [microstrain]')
    ax3.plot(laminate_1b.strains[2][:], t_coo)
    ax3.set_xlabel('e_66 [microstrain]')

    ax4.plot(laminate_1b.stresses[0][:], t_coo)
    ax4.set_ylabel('z (mm)')
    ax4.set_xlabel('s_11 [MPa]')
    ax5.plot(laminate_1b.stresses[1][:], t_coo)
    ax5.set_xlabel('s_22 [MPa]')
    ax6.plot(laminate_1b.stresses[2][:], t_coo)
    ax6.set_xlabel('s_66 [MPa]')

    for ax in axis:
        if len(laminate_1b.lay_up) % 2 == 0:
            print('sym')
            ax.axhline(y = 0, color = 'k', linestyle = '-')
        if len(laminate_1b.lay_up) % 2 != 0:
            print('asym')
            ax.axhline(y = -1*t_coo[int((len(t_coo)-1)/2)], color = 'k', linestyle = '-')
        for t in t_coo:
            ax.axhline(y = t, color = 'k', linestyle = '-')

    plt.show()

if q_2a:
    # initialize material and a lamina
    #material_2a  = Material('UD material', 145.3e3, 8.5e3, 4.58e3, 0.31, None, [1932, 1480, 108, 220, 132.8], 0.2, 230e3) # MPa, [-]
    material_2a  = Material('UD material', 140e3, 10e3, 5e3, 0.3, None, [1500, 1200, 50, 250, 70], 0.2, 230e3) # MPa, [-]
    lamina_2a    = Lamina(material_2a, None, None, None, None, None, None, None, 0)
    #
    lay_up_2a = Layup([0,90], 1, [0.125])
    print(lay_up_2a.stack)
    # initialize laminate
    laminate_2a = Laminate('1b', lay_up_2a.stack, None, np.zeros((3,3) ,dtype='O'), np.zeros((3,3),dtype='O'), np.zeros((3,3),dtype='O'), np.empty((6,6),dtype='O'), np.empty((6,6),dtype='O'), 0, 0, 0, 0, 0, np.empty((3,len(lay_up_2a.stack))), np.empty((3,len(lay_up_2a.stack))))
    # assign lamina orientations and thicknesses based upon specified lay-up 
    GenerateLaminate(lamina_2a, lay_up_2a, laminate_2a)
    # calculate ABD and equivalent properties
    ABD_Matrix(lay_up_2a, laminate_2a)
    ABD_undamaged = laminate_2a.ABD
    # select criteria
    Puck = False
    MaxStress = True
    # initialize force vector F, containing force and moment intensities
    F_2a = np.zeros((6,1)) # (Nx, Ny, Ns, Mx, My, Ms)

    # initialize lists to store failure data
    fpf = [[],[]] # s2, s6
    lpf = [[],[]] # s2, s6
    'FAILURE LOOP'    
    for theta in np.arange(0, 1): #181, 1):
        print('###### UPDATE ANGLE #####')
        print('theta:', theta)
        # assign lamina orientations and thicknesses based upon specified lay-up 
        GenerateLaminate(lamina_2a, lay_up_2a, laminate_2a)
        # calculate ABD and equivalent properties
        ABD_Matrix(lay_up_2a, laminate_2a)
        ABD = laminate_2a.ABD
        #print(ABD)
        if ABD.all() == ABD_undamaged.all():
            # count failures
            failures_occured = 0
            # check failure states
            failure_states = []
            for lamina in laminate_2a.laminae:
                failure_states.append(lamina.failure_state)
            
            if not any(failure_state > 1 for failure_state in failure_states):
                Ny = 100 # N/mm # INITIAL TRANSVERSE LOAD
                while failures_occured != len(lay_up_2a.stack):
                    Ns = Ny*np.tan(np.deg2rad(theta))
                    print()
                    print('###### UPDATE LOAD #####')
                    print('Ny, Ns:', Ny, Ns)
                    F_2a[1], F_2a[2] = Ny, Ns
                    Stress(laminate_2a, F_2a)
                    FI_laminas = []
                    for i, lamina in enumerate(laminate_2a.laminae):
                        print('LAMINA', i+1, ':', lay_up_2a.stack[i], 'DEG')
                        FI  = MaxStressFailure(lamina.stresses, material_2a) # returns highest FI for given lamina
                        print(lamina.stresses)
                        print('FI_lamina:', FI)
                        print()
                        FI_laminas.append(max(FI.values())) # collect all FI
                    print(FI_laminas)
                    FI_max = max(FI_laminas) 
                    # check for failure
                    if FI_max < 0.999:
                        Ny = Ny/FI_max # no failure, so increase load

                    else:
                        print()
                        print('!!!!!!!!!!!!!!!!!!!')
                        print('PLY FAILURE OCCURED')
                        print('!!!!!!!!!!!!!!!!!!!')
                        print()
                        for i, FI in enumerate(FI_laminas):
                            if FI >= 0.999:
                                failures_occured += 1
                                laminate_2a.laminae[i].failure_state = 1
                                # store first and last ply failure
                                if failures_occured == 1 and laminate_2a.laminae[i].failure_state == 1:
                                    fpf[0].append(laminate_2a.laminae[i].stresses[1][0]) # FIX TRANSFORMATION
                                    fpf[1].append(laminate_2a.laminae[i].stresses[2][0])
                                if failures_occured == len(lay_up_2a.stack) and laminate_2a.laminae[i].failure_state == 1:
                                    lpf[0].append(laminate_2a.laminae[i].stresses[1][0])
                                    lpf[1].append(laminate_2a.laminae[i].stresses[2][0])
                                                           
                        if failures_occured != len(lay_up_2a.stack):
                            for i, lamina in enumerate(laminate_2a.laminae):
                                if lamina.failure_state != 0:
                                    'DEGRADATION ~ LEE'
                                    lamina.material.E1  = 0
                                    lamina.material.E2  = 0
                                    lamina.material.G12 = 0
                                    lamina.material.v12 = 0
                                    lamina.material.v21 = 0
                                    laminate_2a.laminae[i] = lamina
                            ABD_Matrix(lay_up_2a, laminate_2a)

                    print('# failures=', failures_occured)

    print(fpf, lpf)            




