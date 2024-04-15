from ABD_Killian import Material, Lamina, Layup, Laminate, GenerateLaminate, ABD_Matrix, T
from StressAnalysis_Killian import Stress
from FailureCriteria_Killian import FibreFailure, MatrixFailure, MaxStressFailure, DegradationRule
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from tqdm import tqdm 
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
    n = 2 # total number of plies = 2*n*5
    t_1a = degradation_factor25 # wasn't given in assignment, assume this value and continue
    # initialize material and a lamina
    material_1a  = Material('UD material', 145.3e3, 8.5e3, 4.58e3, 0.31, None, [1932, 0, 108, 0, 132.8], 0.2, 230e3) # MPa, [-]
    #print(material_1a.v21)
    'TEST TO IMPLEMENT UNCERTAINTIES'
    uncertain = False
    #E1  = ufloat(145.3e3, 3.28e3)
    #E2  = ufloat(8.5e3, 1.28e3)
    #v12 = ufloat(0.31, 0.018)
    #material_1a  = Material('UD material', E1, E2, 4.58e3, v12, None, [1932, 0, 108, 0, 132.8], 0.2, 230e3) # MPa, [-]
    lamina_1a    = Lamina(material_1a, None, None, None, None, None, None, None, 0)

    n = 3
    theta_range = np.arange(0,91,1)
    n_range     = np.arange(1, n+1)
    Ex, Ey, Gxy, vxy, vyx = [],[],[],[],[]
    size = (n_range.size, theta_range.size)
    E1b, E2b, G12b = np.zeros(size), np.zeros(size), np.zeros(size)
    v12b, v21b     = np.zeros(size), np.zeros(size)
    for theta in theta_range:
        for n in n_range:
            base_stack = [15, theta, -theta, 75, 75]
            half_stack = base_stack
            for i in np.arange(1,n):
                half_stack += base_stack 
            lay_up_1a = Layup(half_stack, 1, [t_1a])
            h = lay_up_1a.total_thickness
            # initialize laminate
            laminate_1a = Laminate('1a', lay_up_1a.stack, None, np.zeros((3,3) ,dtype='O'), np.zeros((3,3),dtype='O'), np.zeros((3,3),dtype='O'), np.empty((6,6)), np.empty((6,6)), 0, 0, 0, 0, 0, np.empty((3,len(lay_up_1a.stack))) ,np.empty((3,len(lay_up_1a.stack))))
            # assign lamina orientations and thicknesses based upon specified lay-up 
            GenerateLaminate(lamina_1a, lay_up_1a, laminate_1a) #, laminae)
            # calculate ABD and equivalent properties
            ABD_Matrix(lay_up_1a, laminate_1a)
            # in plane engineering constants
            if n == 1:
                Ex.append(laminate_1a.Ex) 
                Ey.append(laminate_1a.Ey)
                Gxy.append(laminate_1a.Gxy)
                vxy.append(laminate_1a.vxy)
                vyx.append(laminate_1a.vyx)
            # out of plane engineering constants
            d = np.linalg.inv(laminate_1a.D.astype('float64'))
            E1b[n-1][theta]  = (12/(d[0][0]*h**3))
            E2b[n-1][theta]  = (12/(d[1][1]*h**3))
            G12b[n-1][theta] = (12/(d[2][2]*h**3))
            v12b[n-1][theta] = (-d[0][1]/d[0][0])
            v21b[n-1][theta] = (-d[0][1]/d[1][1])

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 5))
    'In Plane Engineering Constants'
    if uncertain: 
        'Upper bound'
        ax1.plot(theta_range, [(i.n + i.s) / 1000 for i in Ex],  color='r')
        ax1.plot(theta_range, [(i.n + i.s) / 1000 for i in Ey],  color='r')
        ax1.plot(theta_range, [(i.n + i.s) / 1000 for i in Gxy], color='r')
        'Mean' 
        ax1.plot(theta_range, [i.n / 1000 for i in Ex],  label = 'Ex')
        ax1.plot(theta_range, [i.n / 1000 for i in Ey],  label = 'Ey')
        ax1.plot(theta_range, [i.n / 1000 for i in Gxy], label = 'Gxy')
        'Lower bound'
        ax1.plot(theta_range, [(i.n - i.s) / 1000 for i in Ex],  color='r')
        ax1.plot(theta_range, [(i.n - i.s) / 1000 for i in Ey],  color='r')
        ax1.plot(theta_range, [(i.n - i.s) / 1000 for i in Gxy], color='r')
    else:
        ax1.plot(theta_range, [i/ 1000 for i in Ex],  label = 'Ex')
        ax1.plot(theta_range, [i/ 1000 for i in Ey],  label = 'Ey')
        ax1.plot(theta_range, [i/ 1000 for i in Gxy], label = 'Gxy')

    'Labelling etc'
    ax1.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_{}'.format(2*n))
    ax1.set_ylabel('GPa')
    ax1.grid()
    ax1.legend()
    if uncertain:
        'Upper bound'
        ax2.plot(theta_range, [i.n + i.s for i in vxy], color='r')
        ax2.plot(theta_range, [i.n + i.s for i in vyx], color='r')
        'Mean'
        ax2.plot(theta_range, [i.n  for i in vxy], label = 'vxy')
        ax2.plot(theta_range, [i.n  for i in vyx], label = 'vyx')
        'Lower bound'
        ax2.plot(theta_range, [i.n - i.s  for i in vxy], color='r')
        ax2.plot(theta_range, [i.n - i.s  for i in vyx], color='r')
    else:
        ax2.plot(theta_range, [i  for i in vxy], label = 'vxy')
        ax2.plot(theta_range, [i  for i in vyx], label = 'vyx')

    ax2.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_{}'.format(2*n))
    ax2.set_ylabel('[-]')
    ax2.grid()
    ax2.legend()

    'Out of Plane Engineering Constants'
    for n in n_range:
        ax3.plot(theta_range, [i / 1000 for i in E1b[n-1]],  label = 'E1b_n={}'.format(n))
        ax3.plot(theta_range, [i / 1000 for i in E2b[n-1]],  label = 'E2b_n={}'.format(n))
        ax3.plot(theta_range, [i / 1000 for i in G12b[n-1]], label = 'G12b_n={}'.format(n))
    ax3.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_{}'.format(2*n))
    ax3.set_ylabel('GPa')
    ax3.grid()
    ax3.legend()
    for n in n_range:
        ax4.plot(theta_range, [i  for i in v12b[n-1]], label = 'v12b_n={}'.format(n))
        ax4.plot(theta_range, [i  for i in v21b[n-1]], label = 'v21b_n={}'.format(n))
    ax4.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_{}'.format(2*n))
    ax4.set_ylabel('[-]')
    ax4.grid()
    ax4.legend()

    plt.show()

'1b' '''STILL NEED TO CONVERT STRESSES TO PRINCIPAL COO SYSTEM: DONE AND VERIFIED WITH CLASS ROOM EXAMPLE'''
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

    ax1.plot(laminate_1b.strains[0][:], t_coo, color = 'r')
    ax1.set_ylabel('z (mm)')
    ax1.set_xlabel('e_11 [microstrain]')
    ax2.plot(laminate_1b.strains[1][:], t_coo, color = 'r')
    ax2.set_xlabel('e_22 [microstrain]')
    ax3.plot(laminate_1b.strains[2][:], t_coo, color = 'r')
    ax3.set_xlabel('e_66 [microstrain]')

    ax4.plot(laminate_1b.stresses[0][:], t_coo, color = 'r')
    ax4.set_ylabel('z (mm)')
    ax4.set_xlabel('s_11 [MPa]')
    ax5.plot(laminate_1b.stresses[1][:], t_coo, color = 'r')
    ax5.set_xlabel('s_22 [MPa]')
    ax6.plot(laminate_1b.stresses[2][:], t_coo, color = 'r')
    ax6.set_xlabel('s_66 [MPa]')
    print(t_coo[0], t_coo[-1])
    for ax in axis:
        ax.vlines(x = 0, color = 'k', linestyle = '--', ymin = -0.3125, ymax = 0.3125)
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
    material_2a  = Material('UD material', 145.3e3, 8.5e3, 4.58e3, 0.31, None, [1932, 1480, 108, 220, 132.8], 0.2, 230e3) # MPa, [-] 
    #material_2a  = Material('UD material', 140e3, 10e3, 5e3, 0.3, None, [1500, 1200, 50, 250, 70], 0.2, 230e3) # MPa, [-]
    lamina_2a    = Lamina(material_2a, None, None, None, None, None, None, None, 0, None)
    # generated a few lay ups for testing:
    lay_up_1  = Layup([0,90, 45, -45, 0, 90, 45, -45], 1, [0.125])
    lay_up_2  = Layup([0], 0, [0.125]) 
    lay_up_exercise = Layup([0, 45, -45, 90],  1, [0.125])

    'ALSO SELECT APPROPRIATE MATERIAL: ASSIGNMENT OR EXERCISE'
    # select the lay up to be used in the analysis:
    lay_up_2a =  lay_up_2
    print(lay_up_2a.stack)
    h = lay_up_2a.total_thickness
    # initialize laminate
    laminate_2a = Laminate('2a', lay_up_2a.stack, None, np.zeros((3,3) ,dtype='O'), np.zeros((3,3),dtype='O'), np.zeros((3,3),dtype='O'), np.zeros((6,6),dtype='O'), np.zeros((6,6),dtype='O'), 0, 0, 0, 0, 0, np.empty((3,len(lay_up_2a.stack))), np.empty((3,len(lay_up_2a.stack))))
    # assign lamina orientations and thicknesses based upon specified lay-up 
    GenerateLaminate(lamina_2a, lay_up_2a, laminate_2a)
    # calculate ABD and equivalent properties
    ABD_Matrix(lay_up_2a, laminate_2a)
    ABD_undamaged = laminate_2a.ABD
    # select criteria
    Puck      = False
    MaxStress = True
    # initialize force vector F, containing force and moment intensities
    F_2a = np.zeros((6,1)) # (Nx, Ny, Ns, Mx, My, Ms)

    print_sanity_check = True
    # initialize lists to store failure data
    fpf = [[],[]] # s2, s6
    lpf = [[],[]] # s2, s6
    NyNs_check = []
    ply_failures = []
    check_mode1 = []
    check_other_modes = []
    'FAILURE LOOP'
    dtheta = 1 
    theta_range = np.arange(58, 59 + dtheta, dtheta)   
    print('yoo')
    for theta in tqdm(theta_range):

        Ny_theta = []
        Ns_theta = []
        sigma_y_theta = []
        tau_xy_theta  = []
        if print_sanity_check:
            print('###### UPDATE ANGLE #####')
            print('theta:', theta)
        # assign lamina orientations and thicknesses based upon specified lay-up 
        GenerateLaminate(lamina_2a, lay_up_2a, laminate_2a)
        # calculate ABD and equivalent properties
        ABD_Matrix(lay_up_2a, laminate_2a)
        ABD = laminate_2a.ABD
        #print(ABD)
        if ABD.all() == ABD_undamaged.all(): # this check is performed to ensure the new angle starts with a fresh laminate
            # count failures
            failures_occured = 0             # failure counter reset 
            # check failure states
            failure_states = []                                # 
            for lamina in laminate_2a.laminae:                 # collecting failure states of all lamina
                failure_states.append(lamina.failure_state)    # 
            
            if all(failure_state == 0 for failure_state in failure_states): # performing a check to make sure no lamina is failed in the fresh stack
                #print('we good to go')
                N = 10 # N/mm # INITIAL TRANSVERSE LOAD
                while failures_occured != len(lay_up_2a.stack):
                    #print('failures_occured:', failures_occured)
                    Ny = N*np.cos(np.deg2rad(theta))
                    Ns = N*np.sin(np.deg2rad(theta))
                    if print_sanity_check:
                        print()
                        print('###### APPLIED LOADS #####')
                        print('Ny, Ns:', Ny, Ns)
                    'APPLY LOADS'
                    F_2a[1], F_2a[2] = Ny, Ns # load applied correctly
                    'STRESS ANALYSIS'
                    Stress(laminate_2a, F_2a) # stresses calculated in principal coo system
                    FI_laminas = []           # list to collect critical F.I.
                    'COLLECT FAILURE INDICATORS'
                    for i, lamina in enumerate(laminate_2a.laminae):
                        if print_sanity_check:
                            print('LAMINA', i+1, ':', lay_up_2a.stack[i], 'DEG')
                        FI  = MaxStressFailure(lamina.stresses, material_2a) # returns dict with all 3 F.I.
                        if print_sanity_check: 
                            #print(lamina.stresses)
                            print('FI_lamina:', FI)
                            print()
                        FI_laminas.append(max(FI.items(), key =lambda x:x[1])) # for every lamina, store the max F.I.
                    if print_sanity_check:    
                        print(FI_laminas)
                    FI_max = max(FI_laminas, key=lambda x:x[1])[1] # take max F.I. accros all lamina to determine load increment 

                    'EVALUATE FAILURE CRITERIA'
                    # check for failure
                    if FI_max < 0.999:   # if max F.I. did not indicate failure
                        if print_sanity_check:
                            print('### LOAD UPDATED ###')
                        N = N/FI_max     # no failure, so increase load
            
 
                    else: # failure occured, so max F.I. = 1 (0.999) 
                        if print_sanity_check:
                            print()
                            print('!!!!!!!!!!!!!!!!!!!')
                            print('PLY FAILURE OCCURED')
                            print('!!!!!!!!!!!!!!!!!!!')
                            print()
                        for i, FI in enumerate(FI_laminas): # loop over maximum F.I. of every lamina
                            lamina = laminate_2a.laminae[i] # take corresponding lamina 
                            if FI[1] >= 0.999: # failure happened 
                                #print('failure registered')      
                                failures_occured += 1
                                lamina.failure_state += 1 # failed 
                                lamina.failure_mode  = FI[0] # store failure mode
                                # transform back to global coo sys
                                T_stress = np.linalg.inv(T(lamina.orientation)[0]) # [1:][:,1:]
                                # store load applied
                                Ny_theta.append(Ny) # storing all Ny at failure for every angle
                                Ns_theta.append(Ns) # idem for Ns

                                # store first and last ply failure
                                if failures_occured == 1 and lamina.failure_state == 1:
                                    sigma_y_theta.append((T_stress @ lamina.stresses)[1][0])
                                    tau_xy_theta.append((T_stress @ lamina.stresses)[2][0])
                                    fpf[0].append((T_stress @ lamina.stresses)[1][0]) # FIX TRANSFORMATION
                                    fpf[1].append((T_stress @ lamina.stresses)[2][0])
                                if failures_occured == len(lay_up_2a.stack) and lamina.failure_state == 1:
                                    sigma_y_theta.append((T_stress @ lamina.stresses)[1][0])
                                    tau_xy_theta.append((T_stress @ lamina.stresses)[2][0])
                                    lpf[0].append((T_stress @ lamina.stresses)[1][0])
                                    lpf[1].append((T_stress @ lamina.stresses)[2][0])
                                                           
                        if failures_occured != len(lay_up_2a.stack):
                            for i, lamina in enumerate(laminate_2a.laminae):
                            # 'IMPLEMENTATION THROUGH FUNCTION'
                            #          if lamina.failure_state == 1:
                            #             lamina = DegradationRule(lamina, lamina.failure_mode, lamina.failure_state)
                            #             #laminate_2a.laminae[i]
                            # ABD_Matrix(lay_up_2a, laminate_2a)
                                    
                                'LEE, straightforward'    
                                if lamina.failure_state == 1: #TODO: fix when degradation rule gets applied (once, multiple times etc)
                                    'DEGRADATION ~ LEE'
                                    degradation_factor = 0  
                                    lamina.material.E1  *= degradation_factor
                                    lamina.material.E2  *= degradation_factor
                                    lamina.material.G12 *= degradation_factor
                                    lamina.material.v12 *= degradation_factor
                                    lamina.material.v21 *= degradation_factor
                                # if lamina.failure_state == 2: #TODO: fix when degradation rule gets applied (once, multiple times etc)
                                #     'DEGRADATION ~ LEE'
                                #     degradation_factor = 0
                                #     lamina.material.E1  *= degradation_factor
                                #     lamina.material.E2  *= degradation_factor
                                #     lamina.material.G12 *= degradation_factor
                                #     lamina.material.v12 *= degradation_factor
                                #     lamina.material.v21 *= degradation_factor
                            ABD_Matrix(lay_up_2a, laminate_2a)
                            print(laminate_2a.ABD)
                            #     'DIRECT IMPLEMENTATION'
                            #     failure_mode = lamina.failure_mode
                            #     failure_state = lamina.failure_state
                            #     value = 2
                            #     #print(failure_mode)
                            #     if failure_state == 1:
                            #         if failure_mode == 'FI1':
                            #             check_mode1.append(failure_mode)
                            #             lamina.material.E1  *= value
                            #             lamina.material.E2  *= value
                            #             lamina.material.G12 *= value
                            #             lamina.material.v12 *= value
                            #             lamina.material.v21 *= value
                            #         if failure_mode != 'FI1':
                            #             check_other_modes.append(failure_mode)
                            #             lamina.material.E2  *= value #0.1
                            #             lamina.material.v21 *= value #0.1
                            #     else:
                            #         lamina.material.E1  *= value
                            #         lamina.material.E2  *= value
                            #         lamina.material.G12 *= value
                            #         lamina.material.v12 *= value
                            #         lamina.material.v21 *= value
                            # ABD_Matrix(lay_up_2a, laminate_2a)
                                
                    #if print_sanity_check:
                    #print('# failures=', failures_occured)

        NyNs_check.append([Ny_theta, Ns_theta])
        ply_failures.append([sigma_y_theta, tau_xy_theta])
    print(check_mode1)
    print(check_other_modes)               
    # control plotting
    plot_stresses = False
    plot_forces   = True
    
    data = NyNs_check
    plt.grid()
    plt.axhline(y=0, color = 'k') #, xmin = -500, xmax = 500)
    plt.axvline(x=0, color = 'k') #, ymin = -500, ymax = 500)
    for theta in np.arange(len(theta_range)):
        if theta == 0:
            plt.scatter(data[theta][0][0], data[theta][1][0],  color = 'b', label = 'first ply failure')
            #plt.scatter(data[theta][0][-1] , data[theta][1][-1] , color = 'r', marker = 'x', label = 'last ply failure')
        else:
            plt.scatter(data[theta][0][0], data[theta][1][0],  color = 'b')
            #plt.scatter(data[theta][0][-1] , data[theta][1][-1] , color = 'r', marker = 'x')

        #plt.plot([0, data[theta][0][-1]], [0, data[theta][1][-1]], linestyle = '--', color = 'k')
        plt.xlabel('N_y [MPa]')
        plt.ylabel('N_s [MPa]')
        plt.legend()

    
    plt.show()




