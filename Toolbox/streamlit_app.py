import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

if columns_to_plot:
    # Create the plot
    plt.figure(figsize=(8, 5))
    for column in columns_to_plot:
        plt.plot(df['X'], df[column], label=column)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Graph of Selected Columns")
    plt.legend()
    plt.grid(True)

    # Display the plot
    st.pyplot(plt)
else:
    st.write("No columns selected for plotting.")
