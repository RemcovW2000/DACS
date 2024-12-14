import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Airfoil import Airfoil
from Laminate import LaminateBuilder
from Member import Member
import copy
import sys
import os

# Add the root directory (DACS) to the Python path
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

#-----------------------------------------------------------------------------------------------------------------------
# streamlit run Toolbox\streamlit_app.py
#-----------------------------------------------------------------------------------------------------------------------

Laminate = LaminateBuilder([45, -45, 0 ,90], True, True, 1)
Member = Member(Laminate)
sparlocations = [80, 250]
sparmembers = [copy.deepcopy(Member) for _ in range(len(sparlocations))]
topmembers = [copy.deepcopy(Member) for _ in range(len(sparlocations)+1)]
botmembers = [copy.deepcopy(Member) for _ in range(len(sparlocations)+1)]

reinforcementpaneltop = LaminateBuilder([45, -45, 0, 0, 0, 0, 0, 0, 90], True, True, 1)
reinforcementpanelbot = LaminateBuilder([45, -45, 0, 0, 0, 0, 0, 0, 90], True, True, 1)
reinforcement_start = 60
reinforcement_end = 100

type = 'NACA2410'
airfoil = Airfoil(type, 1, 300, sparlocations, topmembers, botmembers, sparmembers, reinforcementpaneltop, 60, 100, reinforcementpanelbot, 60, 100)
airfoil.SolveStresses_CSA([30000, 0], [0, 80], [300 / 4, 0])
airfoil.FailureAnalysis_CSA()


# Create sample data
data = {
    'X': np.linspace(0, 10, 50),
    'Y1': np.sin(np.linspace(0, 10, 50)),
    'Y2': np.cos(np.linspace(0, 10, 50))
}
df = pd.DataFrame(data)

# Streamlit App
st.title("Data and Graph Display with Streamlit")

# Display DataFrame
st.header("Sample Data")
st.write("Below is the sample data:")
st.dataframe(df)  # Display as an interactive table

# Plotting
st.header("Graph")
st.write("Select which data to display in the graph:")
columns_to_plot = st.multiselect("Select columns to plot:", options=df.columns[1:], default=['Y1'])
# Create tabs
tab1, tab2 = st.tabs(["Tab 1", "Tab 2"])

# Content for Tab 1
with tab1:
    st.header("Design overview")
    st.write("On this page, an overview of the design is given, with a model airfoil as well as the wing planform.")
    plot = airfoil.plotairfoil()

    # Display the plot in Streamlit
    st.pyplot(plot)

# Content for Tab 2
with tab2:
    st.header("Welcome to Tab 2")
    st.write("This is the content for Tab 2.")
    st.line_chart([1, 2, 3, 4, 5])  # Example chart in Tab 2