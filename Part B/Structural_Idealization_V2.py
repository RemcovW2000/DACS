import numpy as np
import copy
import matplotlib.pyplot as plt
import pandas as pd
from Fuselage import Fuselage
from Boom import Boom 
from Skin import Skin 
from Stringer import Stringer
from Panel import Panel
import scipy.optimize    as opt

#d = 6e3 # [mm] 

np.set_printoptions(linewidth=500, precision=5) #TODO: fix precision

def Structural_Idealization_V2(Mx, Vy, diameter, frame_spacing, stringers, skins, n_custom_stringers):
    n_j = len(skins)     # [-]
    n_s = len(stringers) # [-]

    # print summary
    summary = True
    # specs fuselage section:
    fuselage = Fuselage(diameter, frame_spacing, n_j, n_s - n_custom_stringers)   
    fuselage.stringers = stringers
    fuselage.skins     = skins   

    #NOTE: AFTER THE FIRST LAMINATE DESIGNS, SEPARATE LISTS ARE NEEDED TO ASSIGN SPECIFIC AREAS, MODULI AND THICKNESSES

    'Defining all physical elements'
    theta_j_array  = fuselage.theta_joints_array
    theta_s_array  = fuselage.theta_stringers_array

    # generate all stringers:
    for theta, stringer in zip(theta_s_array, stringers): 
        stringer.location = theta

    # generate all skins:
    for i, skin in enumerate(skins):
        skin.start = theta_j_array[i]
        if i != len(skins) - 1:
            skin.stop = theta_j_array[i + 1]
        if i == len(skins) - 1:
            skin.stop = theta_j_array[0]

        # calculating static moment skin
        skin.CalculateA()
        skin.StaticMomentSkin()

    y_bar                 = fuselage.CalculateNeutralAxis()
    fuselage.neutral_axis = y_bar

    'Defining all ideal elements'
    theta_b_array = np.unique(np.sort(np.concatenate((theta_j_array, theta_s_array))))
    fuselage.theta_booms_array = theta_b_array
    # initialize a boom and a panel:
    boom  = Boom()
    panel = Panel()
    # generate all booms but with zero area:
    for theta in theta_b_array:
        boom = copy.deepcopy(boom)
        boom.location = theta
        fuselage.booms.append(boom)

    # generate all panels with zero thickness:
    for i, theta in enumerate(theta_b_array):
        panel = copy.deepcopy(panel)
        panel.start = theta # first panel starts at first boom onwards
        if theta != theta_b_array[-1]:
            panel.stop  = theta_b_array[i + 1]
        else:
            panel.stop  = theta_b_array[0]
        panel.depth = fuselage.frame_spacing
        fuselage.panels.append(panel)

    # retrieve all panels and booms
    panels    = fuselage.panels
    booms     = fuselage.booms

    # assigning all panels lengths and thicknesses:
    theta_j_array = np.append(np.insert(theta_j_array, 0, 0), 360)

    for panel in panels:
        if panel != panels[-1]:
            panel.length = (np.radians(panel.stop) - np.radians(panel.start))*diameter/2
        else:
            panel.length = (2*np.pi - (np.radians(panel.start) - np.radians(panel.stop)))*diameter/2
        'corresponding to which skin a panel belongs, a thickness and Ex must be assigned'
        for i in range(len(theta_j_array) - 1):
            if theta_j_array[i] <= panel.start < theta_j_array[i + 1]:
                index = i
                if index != 0:
                    panel.Laminate   = skins[index - 1].Laminate
                    panel.thickness  = panel.Laminate.h
                    panel.Ex         = panel.Laminate.Ex
                elif index ==0:
                    panel.Laminate   = skins[- 1].Laminate
                    panel.thickness  = panel.Laminate.h
                    panel.Ex         = panel.Laminate.Ex
                break

        for boom in booms:
            if panel.start == boom.location:
                boom_1 = boom
                y_1    = np.sin(np.radians(boom_1.location))*(diameter/2)
        
            if panel.stop == boom.location:
                boom_2 = boom
                y_2    = np.sin(np.radians(boom_2.location))*(diameter/2)

        panel.B_1 = panel.thickness*panel.length*(2 + (y_2 - y_bar)/(y_1 - y_bar))
        panel.B_2 = panel.thickness*panel.length*(2 + (y_1 - y_bar)/(y_2 - y_bar))

    # assigning all booms an area, depending on contribution of skins and stringers:
    for boom in booms:
        # add area contribution if stringer present at boom location
        if theta_s_array.size > 0:
            if boom.location in theta_s_array:
                stringer_area = stringers[np.where(boom.location == theta_s_array)[0][0]].area
                stringer_Ex   = stringers[np.where(boom.location == theta_s_array)[0][0]].Ex
                boom.area    += stringer_area
                boom.EA      += stringer_Ex*stringer_area
        # add skin area contribution from neighbouring skins
        for panel in panels:
            if panel.start == boom.location:
                panel_1 = panel
            if panel.stop == boom.location:
                panel_2 = panel

        boom.area += panel_1.B_1 + panel_2.B_2
        boom.EA   += panel_1.B_1*panel_1.Ex 
        boom.EA   += panel_2.B_2*panel_2.Ex

    # calculation equivalent properties:
    Ixx = 0
    for boom in booms:
        Ixx += boom.area * (np.sin(np.radians(boom.location))*diameter/2 - y_bar)**2
    # store
    fuselage.Ixx = Ixx

    def AxialStress(Mx, boom, boom_index):
        # boom stress
        boom.sigma = (Mx * np.sin(np.radians(boom.location))*diameter/2 )/(Ixx)
    
        #TODO: NEEDS TO BE FIXED: DONE, VERIFY WITH REMCO
        boom_force = boom.sigma*boom.area
        # distrubiting some stress to stringer
        if theta_s_array.size > 0:
            for stringer in stringers:
                if boom.location == stringer.location:
                    stringer.Fx    = (stringer.area*stringer.Ex/boom.EA)*boom_force
                    stringer.sigma = stringer.Fx/stringer.area

        # distributing stress over panel_1 and panel_2
        panel_1  = panels[boom_index]
        panel_2  = panels[boom_index - 1]

        panel_1.sigma += (panel_1.B_1*panel_1.Ex)*boom_force/(boom.EA)/(panel_1.B_1)
        panel_2.sigma += (panel_2.B_2*panel_2.Ex)*boom_force/(boom.EA)/(panel_2.B_2)

    def ShearFlow(Vy, panel, panel_index):
        for boom in booms:
            if panel.start == boom.location:
                print(Vy, Ixx, np.sin(np.radians(boom.location))*(diameter/2))
                panel.q_b = (-Vy/Ixx)*(boom.area*np.sin(np.radians(boom.location))*(diameter/2)) + panels[panel_index - 1].q_b 
    
    for boom in booms:
        boom.delta_q = (-Vy/Ixx)*(boom.area*np.sin(np.radians(boom.location))*(diameter/2))
        print('boom dq', boom.delta_q)


    panels_to_calculate_q = panels
    #Slice the array to reorder it
    
    for i, panel in enumerate(panels):
        if panel.stop > 90 and panel.start < 90:
            index = i
        if panel.stop == 90:
            index = i
        #break

    shifted = np.concatenate((panels_to_calculate_q[index:], panels_to_calculate_q[:index]))

    print(shifted[0].start)

    for i, panel in enumerate(shifted):
        if i == 0:
            panel.q_b = 0

        else:
            for boom in booms:
                if boom.location == panel.start:
                    panel.q_b = boom.delta_q + shifted[i - 1].q_b
    
        
    for i, boom in enumerate(booms):
        AxialStress(Mx, boom, i)

    # for i, panel in enumerate(panels):
    #     if i == 0:
    #         panel.q_b = 0
    #     else:
    #         ShearFlow(Vy, panel, i)

    # basic shear flows
    #panels_to_calculate_q = panels
    # Slice the array to reorder it
    
    # for i, panel in enumerate(panels):
    #     if panel.stop >= 90:
    #         index = i
    #     break

    # shifted = np.concatenate((panels_to_calculate_q[index:], panels_to_calculate_q[:index]))

    # for i, panel in enumerate(shifted):
    #     if i == 0:
    #         panel.q_b = 0
    #     else:
    #         ShearFlow(Vy, panel, i)




    # constant shear flow to close the cut in panel 1
    moment_q_b = 0
    Am  = np.pi*(diameter/2)**2 # [mm^2]
    for panel in panels:
        moment_q_b += 2*Am*(panel.stop - panel.start)/360

    q_s0 = - moment_q_b/(2*Am) # q_s0 equals zero as shear load is through shear center

    for panel in panels:
        panel.q_b += q_s0 # redundant 

    # sanity check on q performed according to p.653 T.H.G Megson

    # calculating and assigning force intentsities from sigma_z and q
    for panel in panels:
        panel.Nx = panel.sigma*panel.thickness      # TODO: document relation global coo system with local coo system
        panel.Ns = - panel.q_b 

    if summary:
        print('SUMMARY DISCRETIZATION AND LOAD ANALYSIS ----------------------------------------------------------------------------------------')
        print('Skins: ({})'.format(len(skins)) )
        start_skins, stop_skins, t_skins, Ex_skins  = [], [], [], []
        data_skins = {'Start': start_skins, 'Stop': stop_skins, 'Thickness [mm]': t_skins, 'Ex [GPa]': Ex_skins}
        for skin in skins:
            start_skins.append(str(skin.start) + '\u00B0')
            stop_skins.append(str(skin.stop) + '\u00B0')
            t_skins.append(skin.Laminate.h)
            Ex_skins.append(skin.Laminate.Ex/1e3)

        df_skins =pd.DataFrame(data_skins)
        print(df_skins)
        print('---------------------------------------------------------------------------------------------------------------------------------')
        print('Panels: ({}), q_s0 = {}'.format(len(panels), q_s0) )
        start_panels, stop_panels, t_panels, length_panels, sigma_panels, Ex_panels = [],[],[],[],[],[]
        q_panels, Nx_panels, Ns_panels, N_panels = [],[],[],[]
        data_panels = {'Start':start_panels, 'Stop': stop_panels, 'Thickness [mm]': t_panels, 
                    'Length [mm]': length_panels, 'Ex [GPa]': Ex_panels, 'Axial Stress [MPa]': sigma_panels, 
                        'q [N/mm]': q_panels, 'Nx [N/mm]': Nx_panels, 'Ns [N/mm]': Ns_panels,'N [N/mm]': N_panels}
        #TODO: REMOVE EQUIVALENT AREA 
        for panel in panels:
            start_panels.append(str(panel.start) + '\u00B0')
            stop_panels.append(str(panel.stop) + '\u00B0')
            t_panels.append(panel.thickness)
            length_panels.append(round(panel.length,3))
            sigma_panels.append(round(panel.sigma, 3))
            #eq_area.append(round(panel.eq_area, 3))
            q_panels.append(round(panel.q_b - q_s0, 3))
            Nx_panels.append(round(panel.Nx, 3))
            Ns_panels.append(round(panel.Ns, 3))
            Ex_panels.append(panel.Laminate.Ex/1e3)
            N_panels.append(round(np.sqrt(panel.Nx**2 + panel.Ns**2), 3))
            

        df_panels = pd.DataFrame(data_panels)
        print(df_panels)
        if theta_s_array.size > 0:
            print('---------------------------------------------------------------------------------------------------------------------------------')
            print('Stringers: ({})'.format(len(stringers)) )
            location_stringers, area_stringers, sigma_stringers, Ex_stringers = [],[],[],[] #TODO: CHECK UNITS STRESS: N, mm -> MPa
            data_stringers = {'Location': location_stringers, 'Area [mm^2]': area_stringers,
                            'Ex [GPa]':Ex_stringers, 'Axial Stress [MPa]':sigma_stringers, 'Nx [MPa]': sigma_stringers}
            for stringer in stringers:
                location_stringers.append(str(stringer.location) + '\u00B0')
                area_stringers.append(round(stringer.area,3))
                sigma_stringers.append(round(stringer.sigma,3))
                Ex_stringers.append(stringer.Ex/1e3)
            
            df_stringers = pd.DataFrame(data_stringers)
            print(df_stringers)
        print('---------------------------------------------------------------------------------------------------------------------------------')
        print('Booms: ({})'.format(len(booms)) )
        location_booms, area_booms, sigma_booms, EA_booms = [],[], [], []
        data_booms = {'Location:':location_booms, 'Area [mm^2]': area_booms, 'EA [kN]':EA_booms, 'Axial Stress [MPa]': sigma_booms}
        for boom in booms:
            location_booms.append(str(boom.location) + '\u00B0')
            area_booms.append(round(boom.area, 3))
            sigma_booms.append(np.round(boom.sigma,3))
            EA_booms.append(boom.EA/1e3 )

        df_booms = pd.DataFrame(data_booms)
        print(df_booms)

    return fuselage 

