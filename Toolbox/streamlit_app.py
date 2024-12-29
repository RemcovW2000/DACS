import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from airfoil import Airfoil
from laminate import laminate_builder
from member import Member
from wing import generate_chord_and_leading_edge_elyptical, generate_chord_and_leading_edge_tapered,  generate_lift_elyptical, Wing
from Data.Panels import Sandwiches, Laminates
import copy
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Add the root directory (DACS) to the Python path
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

# ----------------------------------------------------------------------------------------------------------------------
# streamlit run Toolbox\streamlit_app.py
# ----------------------------------------------------------------------------------------------------------------------

# Example usage:
n = 200  # Number of points
half_span = 1500  # Half span of the wing
chord_at_root = 350  # Chord length at the root
taper_ratio = 0.4
# chord_lengths, leading_edge_locations = generate_chord_and_leading_edge_elyptical(n, half_span, chord_at_root)
chord_lengths, leading_edge_locations = generate_chord_and_leading_edge_tapered(n, half_span, chord_at_root, taper_ratio)

# Example usage:
weight = 25         # kg
load_factor = 4     # g
g = 9.81            # m/s
lift = weight*load_factor*g
lift_distribution = generate_lift_elyptical(n, half_span, lift)

thickness_distribution = [10, 12, 14, 16, 18]  # Example thickness distribution
half_span = 1500  # Half span of the wing
spar_coordinates = [[[300/4, 0], [15 + 250/4, 1000], [75, 1500]],
                   [[200, 0], [185, 1000], [75, 1500]]]  # Example spar coordinates
rib_coordinates = [0, 600, 1200, 1400]

reinforcementcoordinates = [[300/4, 0], [75, 500], [75, 1000]] # list containing reinforcement coordinates: [[x, y], [x, y] ... ]

wing = Wing(lift_distribution, chord_lengths, leading_edge_locations, thickness_distribution, half_span, rib_coordinates, spar_coordinates)

wing.toppanels = [[Laminates['CFQIezcomposites_spreadtow'], 500], [Laminates['CFQIezcomposites_spreadtow'], 1600]]
wing.botpanels = [[Laminates['CFQIezcomposites_spreadtow'], 500], [Laminates['CFQIezcomposites_spreadtow'], 1600]]
wing.trpanels = [[Laminates['Reinforcement'], 700], [Laminates['Reinforcement'], 1000]]
wing.brpanels = [[Laminates['Reinforcement'], 700], [Laminates['Reinforcement'], 1000]]
wing.trwidth = 50
wing.brwidth = 50
wing.trcoordinates = reinforcementcoordinates
wing.brcoordinates = reinforcementcoordinates


wing.sparpanels = Laminates['CFQIezcomposites_spreadtow']
moment_distribution = wing.internal_moment()
shear_distribution = wing.shear_force()
wing.tip_buffer = 50
wing.nr_airfoils = 5

wing.generate_airfoil_objects()
wing.solve_stresses_CSA()
wing.failure_analysis()
wing.calculate_weight()

# Create sample data
data = {
    'X': np.linspace(0, 10, 50),
    'Y1': np.sin(np.linspace(0, 10, 50)),
    'Y2': np.cos(np.linspace(0, 10, 50))
}
df = pd.DataFrame(data)
# Create tabs
tab1, tab2, tab3 = st.tabs(["Wing design overview", "Airfoil force intensities", "Failure analysis data overview"])

# Content for Tab 1
with tab1:
    st.header("Design overview")
    st.write("On this page, an overview of the design is given, with a model airfoil as well as the wing planform.")
    root_airfoil = wing.airfoils[0]
    tip_airfoil = wing.airfoils[-1]

    st.write(" overview of the wing planform:")
    plotplanform = wing.plot_wing_planform()
    st.pyplot(plotplanform)

    # for middle airfoil, find middle object in list:
    middle_index = (len(wing.airfoils) - 1) // 2
    mid_airfoil = wing.airfoils[middle_index]

    st.write("Root airfoil cross section design:")
    plotroot = root_airfoil.plot_airfoil()
    st.pyplot(plotroot)

    st.write("middle airfoil cross section design:")
    plotmid = mid_airfoil.plot_airfoil()
    st.pyplot(plotmid)

    st.write("tip airfoil cross section design:")
    plottip = tip_airfoil.plot_airfoil()
    st.pyplot(plottip)

# Content for Tab 2
with tab2:
    st.header("Overview of force intensities (N/mm) of cross sections")
    st.write("")
    # Slider for selecting the location
    span_values = np.linspace(0, wing.half_span, len(wing.airfoils))
    slider_value = st.slider("Select a spanwise location:", min_value=0.0, max_value=float(wing.half_span), step=10.0)

    # Find the airfoil closest to the selected location
    airfoilys = [airfoil.y for airfoil in wing.airfoils]
    closest_index = np.argmin(np.abs(np.array(airfoilys) - slider_value))
    closest_airfoil = wing.airfoils[closest_index]

    # Generate and display the plot for the closest airfoil
    st.write(
        f"Displaying airfoil for spanwise location: {slider_value:.2f} mm (Closest to {airfoilys[closest_index]:.2f} mm)")

    st.pyplot(closest_airfoil.plot_airfoil())
    st.pyplot(closest_airfoil.plot_shear_flow())

    # Create DataFrame for display
    Sxgiven = np.round(closest_airfoil.Sx, 2)
    Sygiven = np.round(closest_airfoil.Sy, 2)
    Mgiven = np.round(closest_airfoil.Moment_applied, 2)

    Sxcalc, Sycalc, Mcalc = np.round(closest_airfoil.shear_force_analysis(), 2)
    errors = [abs(Sxcalc - Sxgiven), abs(Sycalc - Sygiven), abs(Mcalc-Mgiven)]
    errorspercent = [abs(Sxcalc - Sxgiven)/(Sxgiven + 1e-9), abs(Sycalc - Sygiven)/(Sygiven+1e-9), abs(Mcalc-Mgiven)/(Mgiven+1e-9)]
    errorspercent = [errorspercent[i] * 100 if errorspercent[i] < 1 else 'high' for i in range(len(errorspercent))]

    row_names = ['x direction shear force (N)',
                 'y direction shear force (N)',
                 'shear moment (Nmm)']
    shear_data = pd.DataFrame({
        "Calculated": [Sxcalc, Sycalc, Mcalc],
        "Given": [Sxgiven, Sygiven, Mgiven],
        "Error": errors,
        "Error (percent)" : errorspercent
    },
    index=row_names
    )

    st.subheader("Shear Flow Comparison Table")
    st.write("Below is a comparison of calculated and given shear flows, along with the error.")
    st.table(shear_data)  # Use st.dataframe if interactivity is needed
    st.pyplot(closest_airfoil.plot_normalforce_intensity())

with tab3:
    st.header("Global performance analysis")
    st.write("This page gives a performance analysis of the entire wing")

    st.write("Total weight: ", wing.weight)

    st.header("Max FI plotted as function of half span:")
    st.pyplot(wing.plot_max_FI())

    st.header("Full overview of failures:")
    st.write(wing.get_hierarchy())