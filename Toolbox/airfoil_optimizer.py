import numpy as np
import numpy as np
import os
import subprocess
from scipy.optimize import minimize
import time


def airfoil_analysis(thickness, thickest_point, camber, Cl, Re, airfoil_name, airfoil_filename, n_iter=1000):
    # Prepare XFOIL input file
    if os.path.exists("polar_file.txt"):
        os.remove("polar_file.txt")

    input_file = open("input_file.in", 'w')
    input_file.write("\n")
    input_file.write("plop \n")
    input_file.write("G \n")
    input_file.write("\n")
    input_file.write("\n")

    input_file.write("naca 2410")
    input_file.write("\n")
    # input_file.write("load " + airfoil_filename + '\n')
    # input_file.write(airfoil_name + '\n')
    # input_file.write("panel\n")

    # Modify the camber, thickness, and thickest point
    input_file.write("GDES\n")
    input_file.write("TSET\n")  # Set thickness factor
    input_file.write(f"{camber}\n")
    input_file.write("\n")
    input_file.write(f"{thickness}\n")
    input_file.write("high\n")  # Modify thickest point location
    input_file.write(f"{thickest_point}\n")
    input_file.write("\n")
    input_file.write("x\n")
    input_file.write("\n")
    input_file.write('mdes \n')
    input_file.write('filt \n')
    input_file.write('filt \n')
    input_file.write('filt \n')
    input_file.write('filt \n')
    input_file.write('filt \n')
    input_file.write('x \n')
    input_file.write('\n')
    input_file.write('pane \n')

    # Aerodynamic analysis
    input_file.write("OPER\n")
    input_file.write(f"v {Re}\n")
    input_file.write("PACC\n")
    input_file.write("polar_file.txt\n\n")
    input_file.write(f"ITER {n_iter}\n")
    input_file.write(f"CSeq {Cl} {Cl} 1\n")
    input_file.write("\n\n")
    input_file.write("quit\n")
    input_file.close()

    # Running XFOIL
    try:
        with open("input_file.in", 'r') as input_file:
            input_contents = input_file.read()

        xfoil_process = subprocess.Popen(["xfoil.exe"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE, text=True)
        xfoil_output, xfoil_error = xfoil_process.communicate(input=input_contents, timeout=2)  # Timeout in seconds
        if xfoil_process.returncode == 0:
            pass
        else:
            print("XFOIL analysis did not converge or encountered an error.")
            print("XFOIL Output:")
            print(xfoil_output)
            print("XFOIL Error:")
            print(xfoil_error)
    except subprocess.TimeoutExpired:
        xfoil_process.terminate()
        print("XFOIL analysis timed out.")

    # Reading results from polar file
    try:
        polar_data = np.loadtxt("polar_file.txt", skiprows=12)
    except ValueError:
        polar_data = np.zeros(3)
        polar_data.fill(2)  # Safety feature in case XFOIL output is empty

    if polar_data.any() == None or len(polar_data) < 3:
        Cd = None
        AOA = None
    else:
        Cd = polar_data[2]
        AOA = polar_data[0]

    return Cd, AOA


# Example usage
Re = 1e6  # Reynolds number
Cl = 0  # Desired lift coefficient
thickness = 0.1  # Thickness as a fraction of chord (e.g., 0.12 for 12%)
thickest_point = 0.34  # Position of thickest point along chord (e.g., 0.3 for 30% chord)
camber = 0.1  # Camber as a fraction of chord (e.g., 0.04 for 4%)
airfoil_name = "e395-il"  # Example airfoil name
airfoil_filename = "e395-il.txt"  # Example airfoil file

start_time = time.time()
Cd, AOA = airfoil_analysis(thickness, thickest_point, camber, Cl, Re, airfoil_name, airfoil_filename)
end_time = time.time()
print(end_time - start_time)
print()
print(f"Minimum Cd: {Cd}")
print(f"Angle of Attack for desired Cl: {AOA}")