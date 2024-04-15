from ABD_Killian import Material, Lamina, Layup, Laminate, GenerateLaminate, ABD_Matrix
from StressAnalysis_Killian import Stress
from FailureCriteria_Killian import MaxStressFailure, PuckFailure, DegradationRule
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from tqdm import tqdm 

# print settings
np.set_printoptions(linewidth=300, precision = 3)

# control what questions are active
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
    t_1a = 0.125 # wasn't given in assignment, assume this value and continue
    # initialize material and a lamina
    material_1a  = Material('UD material', 145.3e3, 8.5e3, 4.58e3, 0.31, None, [1932, 0, 108, 0, 132.8], 0.2, 230e3) # MPa, [-]
    #print(material_1a.v21)
    'TEST TO IMPLEMENT UNCERTAINTIES'
    uncertain = False
    lamina_1a    = Lamina(material_1a, None, None, None, None, None, None, None, 0, None)

    n = 8
    theta_range = np.arange(0,91,1)
    n_range     = np.arange(1, n+1)
    Ex, Ey, Gxy, vxy, vyx = [],[],[],[],[]
    E1b_converged = []
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
            #print(lay_up_1a.stack)
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
            v12b[n-1][theta] = (-d[0][1]/d[1][1])
            v21b[n-1][theta] = (-d[0][1]/d[0][0])

            if theta==0: 
                Q = laminate_1a.laminae[0].Q
                E1b_con = 4*(Q[0][0] + Q[0][1]**2/Q[1][1])
                E1b_converged.append(E1b_con)
        
    plt.plot(n_range, E1b_converged)
    plt.show()

    fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(10, 5)) # , (ax3, ax4)
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
    ax1.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_s')  #_{}'.format(2*n))
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

    ax2.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_s') #) _{}'.format(2*n))
    ax2.set_ylabel('[-]')
    ax2.grid()
    ax2.legend()
    plt.show()

    # 'Out of Plane Engineering Constants'

    fig, ((ax3, ax4, ax5)) = plt.subplots(1, 3, figsize=(10, 5))
    # E1b
    for n in n_range:
        ax3.plot(theta_range, [i / 1000 for i in E1b[n-1]],  label = 'E1b_n={}'.format(n))
    ax3.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_ns')
    ax3.set_ylabel('GPa')
    ax3.grid()
    ax3.legend()
    # E2b
    for n in n_range:
        ax4.plot(theta_range, [i / 1000 for i in E2b[n-1]],  label = 'E2b_n={}'.format(n))
    ax4.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_ns')
    ax4.set_ylabel('GPa')
    ax4.grid()
    ax4.legend()
    # G12b
    for n in n_range:
        ax5.plot(theta_range, [i / 1000 for i in G12b[n-1]], label = 'G12b_n={}'.format(n))
    ax5.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_ns')
    ax5.set_ylabel('GPa')
    ax5.grid()
    ax5.legend()

    fig, ((ax6, ax7)) = plt.subplots(1, 2, figsize=(10, 5))
    for n in n_range:
        ax6.plot(theta_range, [i  for i in v12b[n-1]], label = 'v12b_n={}'.format(n))
    ax6.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_2n')
    ax6.set_ylabel('[-]')
    ax6.grid()
    ax6.legend()

    for n in n_range:
        ax7.plot(theta_range, [i  for i in v21b[n-1]], label = 'v21b_n={}'.format(n))
    ax7.set_xlabel('\u03B8 (deg) from [15, \u03B8, -\u03B8, 75\u2082]_2n')
    ax7.set_ylabel('[-]')
    ax7.grid()
    ax7.legend()


    plt.show()

'1b' '''STILL NEED TO CONVERT STRESSES TO PRINCIPAL COO SYSTEM: DONE AND VERIFIED WITH CLASS ROOM EXAMPLE'''
if q_1b:
    # initialize material and a lamina
    material_1b  = Material('UD material', 145.3e3, 8.5e3, 4.58e3, 0.31, None, [1932, 0, 108, 0, 132.8], 0.2, 230e3) # MPa, [-]
    #material_1b  = Material('UD material', 140e3, 10e3, 5e3, 0.3, None, [1932, 0, 108, 0, 132.8], 0.2, 230e3) # MPa, [-]
    lamina_1b    = Lamina(material_1b, None, None, None, None, None, None, None, 0, None)
    #
    t_1b = [0.125]
    lay_up_1b = Layup([0,0,90,30,90], 0, t_1b) # Layup([0,90], 1, t_1b) #Layup([0,90], 1, t_1b) #
    

    # initialize laminate
    laminate_1b = Laminate('1b', lay_up_1b.stack, None, np.zeros((3,3)), np.zeros((3,3)), 
                           np.zeros((3,3)), np.empty((6,6)), np.empty((6,6)), 0, 0, 0, 0, 0, np.empty((3,len(lay_up_1b.stack))), 
                           np.empty((3,len(lay_up_1b.stack))), np.empty((3,len(lay_up_1b.stack))), np.empty((3,len(lay_up_1b.stack))))
    
    # assign lamina orientations and thicknesses based upon specified lay-up 
    GenerateLaminate(lamina_1b, lay_up_1b, laminate_1b)
    # calculate ABD and equivalent properties
    ABD_Matrix(lay_up_1b, laminate_1b)
    # initialize force vector F, containing force and moment intensities
    F_1b = np.zeros((6,1)) # (Nx, Ny, Ns, Mx, My, Ms)
    # apply loadcase
    F_1b[0]  = 0.2e-1 # N/mm -> results in stresses in MPa (N/mm^2)
    F_1b[1]  = 1.8e1
    F_1b[3]  = 18e3

    #F_1b[0] = 100 #N/mm
   
    Stress(laminate_1b, F_1b)

    t_coo = []
    for lamina in laminate_1b.laminae:
        t_coo.append(max([lamina.z0, lamina.z1], key=abs))

    print(laminate_1b.maxstrains)
    print(laminate_1b.minstrains)
    print()
    print(laminate_1b.maxstresses)
    print(laminate_1b.minstresses)

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3) #, figsize=(10, 5))

    axis = [ax1, ax2, ax3, ax4, ax5, ax6]

    # ax1 -> e11 -> max/min strains[0]
    # ax2 -> e22 -> max/min strains[1]
    # ax2 -> e66 -> max/min strains[2]
    # idem for stresses, ax3, ax4, ax5
    laminae = laminate_1b.laminae
    maxstrains = laminate_1b.maxstrains*1e3
    minstrains = laminate_1b.minstrains*1e3
    maxstresses = laminate_1b.maxstresses/1e3
    minstresses = laminate_1b.minstresses/1e3

    for i in range(len(lay_up_1b.stack)):
        for j, ax in enumerate([ax1, ax2, ax3]):
            ax.plot([maxstrains[j][i], minstrains[j][i]], [laminae[i].z0, (laminae[i].z0 + t_1b)[0]], color = 'r', zorder = 10)
            if i != len(lay_up_1b.stack) - 1:
                ax.plot([minstrains[j][i], maxstrains[j][i+1]], [(laminae[i].z0 + t_1b)[0], laminae[i+1].z0], color = 'r',  zorder = 10)
            if i == 0: 
                ax.plot([0, maxstrains[j][i]], [laminae[i].z0, laminae[i].z0], color = 'r',  zorder = 10)
            if i == len(lay_up_1b.stack) - 1:
                ax.plot([minstrains[j][i], 0], [(laminae[i].z0 + t_1b)[0], (laminae[i].z0 + t_1b)[0]], color = 'r',  zorder = 10)
    
    for i in range(len(lay_up_1b.stack)):
        for j, ax in enumerate([ax4, ax5, ax6]):
            ax.plot([maxstresses[j][i], minstresses[j][i]], [laminae[i].z0, (laminae[i].z0 + t_1b)[0]], color = 'r',  zorder = 10)
            if i != len(lay_up_1b.stack) - 1:
                ax.plot([minstresses[j][i], maxstresses[j][i+1]], [(laminae[i].z0 + t_1b)[0], laminae[i+1].z0], color = 'r',  zorder = 10)
            if i == 0: 
                ax.plot([0, maxstresses[j][i]], [laminae[i].z0, laminae[i].z0], color = 'r',  zorder = 10)
            if i == len(lay_up_1b.stack) - 1:
                ax.plot([minstresses[j][i], 0], [(laminae[i].z0 + t_1b)[0], (laminae[i].z0 + t_1b)[0]], color = 'r',  zorder = 10)

    ax1.set_ylabel('z (mm)')
    ax1.text()
    ax1.set_xlabel('\u03B5_11 [microstrain]')
    ax2.set_xlabel('\u03B5_22 [microstrain]')
    ax3.set_xlabel('\u03B5_66 [microstrain]')
    ax4.set_ylabel('z (mm)')
    ax4.set_xlabel('\u03C3_11 [GPa]')
    ax5.set_xlabel('\u03C3_22 [GPa]')
    ax6.set_xlabel('\u03C3_66 [GPa]')

    print(t_coo[0], t_coo[-1])
    for ax in axis:
        ax.vlines(x = 0, color = 'k', linestyle = '--', ymin = -0.3125, ymax = 0.3125)
        if len(laminate_1b.lay_up) % 2 == 0:
            #print('sym')
            ax.axhline(y = 0, color = 'k', linestyle = '-')
        if len(laminate_1b.lay_up) % 2 != 0:
            #print('asym')
            ax.axhline(y = -1*t_coo[int((len(t_coo)-1)/2)], color = 'k', linestyle = '-')
        for t in t_coo:
            ax.axhline(y = t, color = 'k', linestyle = '-')

    plt.show()

if q_2a:
    # initialize material and a lamina
    material_2a  = Material('UD material', 145.3e3, 8.5e3, 4.58e3, 0.31, None, [1932, 1480, 108, 220, 132.8] , 0.2, 230e3) # MPa, [-]  
    lamina_2a    = Lamina(material_2a, None, None, None, None, None, None, None, 0, None)
    # lay up for assignment:
    lay_up_assignment  = Layup([0,90, 45, -45, 0, 90, 45, -45], 1, [0.125])
    # generated a few lay ups for testing:
    lay_up_2        = Layup([45, 45], 0, [0.125]) 
    lay_up_exercise = Layup([0, 45, -45, 90],  1, [0.125])
    lay_up_test     = Layup([0,90,45,-45,], 0, [0.125]) 
    'USE THIS LAYUP ^^^ with theta 10-11 deg to illustrate concept of damage tolerance'
    lay_up_UD    = Layup([0], 0, [0.125])
    lay_up_Remco = Layup([0,90], 1, [0.125])

    'ALSO SELECT APPROPRIATE MATERIAL: ASSIGNMENT OR EXERCISE'
    # select the lay up to be used in the analysis:
    lay_up_2a =  lay_up_assignment
    print('Laminate:', lay_up_2a.stack)
    print('-----------------------------')
    print('FAILURE ENVELOPE CALCULATION:')
    print('-----------------------------')
    h = lay_up_2a.total_thickness
    # initialize laminate
    laminate_2a = Laminate('2a', lay_up_2a.stack, None, np.zeros((3,3)), np.zeros((3,3)), 
                           np.zeros((3,3)), np.zeros((6,6)), np.zeros((6,6)), 0, 0, 0, 0, 
                           0, np.empty((3,len(lay_up_2a.stack))), np.empty((3,len(lay_up_2a.stack)))
                           , np.empty((3,len(lay_up_2a.stack))), np.empty((3,len(lay_up_2a.stack))))
    # assign lamina orientations and thicknesses based upon specified lay-up 
    GenerateLaminate(lamina_2a, lay_up_2a, laminate_2a)
    # calculate ABD and equivalent properties
    ABD_Matrix(lay_up_2a, laminate_2a)
    ABD_undamaged = laminate_2a.ABD

    'CONTROL PRINTING OF INTERMEDIATE INFORMATION BELOW'
    print_sanity_check = False

    if print_sanity_check:
            print('ABD_undamaged:')
            print(ABD_undamaged)
    # select criteria
    Puck      = True
    MaxStress = False
    strain_Puck = []
    NyNs_Puck   = []
    strain_Max  = []
    NyNs_Max    = []
    failure_criteria = [Puck, MaxStress]
    for criteria in failure_criteria:
        # initialize force vector F, containing force and moment intensities
        F_2a = np.zeros((6,1)) # (Nx, Ny, Ns, Mx, My, Ms)
        # initialize lists to store failure data
        'FAILURE LOOP'
        dtheta = 1 # step in angle theta
        theta_range = np.arange(0, 360 + dtheta, dtheta)

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
            if ABD.all() == ABD_undamaged.all(): # this check is performed to ensure the new angle starts with a fresh laminate
                # count failures
                broken_plies = 0             # failure counter reset 
                # check failure states
                failure_states = []                                  # 
                for lamina in laminate_2a.laminae:                   # collecting failure states of all lamina
                    failure_states.append(lamina.failure_counter)    # 
                
                if all(failure_counter == 0 for failure_counter in failure_states): # performing a check to make sure no lamina is failed in the fresh stack
                    #print('we good to go')
                    N = 10 # N/mm # INITIAL TRANSVERSE LOAD
                    while   broken_plies < len(lay_up_2a.stack):
                        if print_sanity_check:
                            print('# broken_plies =', broken_plies)
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
                            if criteria == MaxStress:
                                FI  = MaxStressFailure(lamina.stresses[0], material_2a) # returns dict with all 3 F.I.
                            if criteria == Puck:
                                FI  = PuckFailure(lamina.stresses[0], material_2a, print_sanity_check) # returns dict with all 4 F.I.
                            if print_sanity_check: 
                                print('FI_lamina:', FI)
                                print()
                            FI_laminas.append(max(FI.items(), key = lambda x:x[1])) # for every lamina, store the max F.I.
                        if print_sanity_check:    
                            print(FI_laminas)
                        FI_max = max(FI_laminas, key = lambda x:x[1]) # take max F.I. accros all lamina to determine load increment 
                        'EVALUATE FAILURE CRITERIA'
                        # check for failure
                        if FI_max[1] < 0.999:   # if max F.I. did not indicate failure
                            if print_sanity_check:
                                print()
                                print('### LOAD UPDATED ###')
                            if MaxStress:
                                N = N/(FI_max[1])     # no failure, so increase load
                            if Puck:
                                N = N/(FI_max[1]) 
                
                        else: # failure occured, so max F.I. = 1 (0.999) 
                            if print_sanity_check:
                                print()
                                print('!!!!!!!!!!!!!!!!!!!')
                                print('PLY DAMAGE OCCURED')
                                print('!!!!!!!!!!!!!!!!!!!')
                                print()
                            for i, FI in enumerate(FI_laminas): # loop over maximum F.I. of every lamina
                                lamina = laminate_2a.laminae[i] # take corresponding lamina 
                                if FI[1] >= 0.999: # failure happened      
                                    lamina.failure_counter += 1 # tracking amount of local failures
                                    lamina.failure_mode     = FI[0] # store failure mode
                                    if lamina.failure_counter == 1 and lamina.failure_mode == 'FI1':
                                        if print_sanity_check:
                                            print('fibre failure registered, LAMINA:', i + 1)
                                        broken_plies       += 1 # tracking amount of global failures
                                    if lamina.failure_counter >= 2:
                                        if print_sanity_check:
                                            print('2nd failure registered, LAMINA:', i + 1)
                                        broken_plies       += 1 # tracking amount of global failures
                                    if (lamina.failure_mode == 'FI2' or lamina.failure_mode == 'FI3') and lamina.failure_counter ==1:
                                        if print_sanity_check:
                                            print('matrix failure registered, LAMINA:', i + 1)

                                    # store load applied
                                    Ny_theta.append(Ny) # storing all Ny at failure for every angle
                                    Ns_theta.append(Ns) # idem for Ns
                                    #strain_Puck.append()

                            if  broken_plies < len(lay_up_2a.stack):
                                for i, lamina in enumerate(laminate_2a.laminae):
                                    if lamina.failure_counter <= 2:
                                        DegradationRule(lamina)
                                ABD_Matrix(lay_up_2a, laminate_2a)
                                if print_sanity_check:
                                    print()
                                    print('Degraded ABD:')
                                    print(laminate_2a.ABD)
                            
                                    
                        if print_sanity_check:
                            print()
                            print('# broken_plies =', broken_plies)
                            print('---------------------------------------------------------------------')
            if criteria == MaxStress:
                NyNs_Max.append([Ny_theta, Ns_theta])
            if criteria == Puck:
                NyNs_Puck.append([Ny_theta, Ns_theta])
             
 

    #data = NyNs #[i/h for i in NyNs]
    plt.grid() 
    plt.axhline(y=0, color = 'k')
    plt.axvline(x=0, color = 'k')

    to_be_plotted = [NyNs_Puck] #NyNs_Max, 

    for i in to_be_plotted:
        data = i
        if i == NyNs_Max:
            color1 = 'g'
            color2 = 'k'
            method = 'Max Stress'
        else:
            color1 = 'b'
            color2 = 'r'
            method = 'Puck'

        for theta in np.arange(len(theta_range)):
            if theta == 0:
                plt.scatter(data[theta][0][0], data[theta][1][0],  color = color1, label = 'first ply DAMAGE with ' + method)
                plt.scatter(data[theta][0][-1] , data[theta][1][-1] , color = color2, marker = 'x', label = 'last ply FAILURE with ' + method)
            else:
                plt.scatter(data[theta][0][0], data[theta][1][0],  color = color1)
                plt.scatter(data[theta][0][-1] , data[theta][1][-1] , color = color2, marker = 'x')

        plt.title('Failure Envelope: Puck Criteria') 
        plt.xlabel('N_y [N/mm]')
        plt.ylabel('N_s [N/mm]')
        plt.legend()
    plt.axis('equal')
    plt.show()

