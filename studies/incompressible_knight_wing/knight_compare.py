import numpy as np
import matplotlib.pyplot as plt
from dev.helper_scripts.paraview_functions import *


if __name__=="__main__":

    # Locations of pressure slices
    # Station 1 : 14.75 in
    # Station 2 : 14.25 in
    # Station 3 : 13.5 in
    # Station 4 : 12.5 in
    # Station 5 : 11.0 in
    # Station 6 : 9.0 in
    # Station 7 : 6.0 in
    # Station 8 : 3.0 in

    # Read in experimental data
    raw_data = np.genfromtxt("studies/incompressible_knight_wing/experimental pressures -8.csv", delimiter=',', skip_header=2)

    # Get pressure corrections
    C_p_refs = raw_data[0,33::2]

    # Clean up data
    C_p_lower = []
    x_lower = []
    C_p_upper = []
    x_upper = []
    for i in range(8):

        # Get pressure correction
        #C_p_ref = raw_data[0,2*i+31]

        # Get data
        x_upper.append(raw_data[:,4*i])
        C_p_upper.append(raw_data[:,4*i+1])
        x_lower.append(raw_data[:,4*i+2])
        C_p_lower.append(raw_data[:,4*i+3])

        # Remove nans and correct pressure
        x_upper[i] = x_upper[i][np.logical_not(np.isnan(x_upper[i]))]
        x_lower[i] = x_lower[i][np.logical_not(np.isnan(x_lower[i]))]
        C_p_upper[i] = C_p_upper[i][np.logical_not(np.isnan(C_p_upper[i]))] - C_p_refs[i]
        C_p_lower[i] = C_p_lower[i][np.logical_not(np.isnan(C_p_lower[i]))] - C_p_refs[i]

        # Load panel data
        panel_data = np.genfromtxt("station {0} morino.csv".format(i+1), delimiter=',', skip_header=1)
        C_p_panel = panel_data[:,0]
        x_panel = panel_data[:,1]
        x_min = min(x_panel)
        x_max = max(x_panel)
        x_panel /= x_max-x_min # Normalize
        x_panel *= -1.0 # Flip
        x_panel += 0.25 # Shift

        # plot
        plt.figure(figsize=(5.5, 5))
        plt.plot(x_upper[i], C_p_upper[i], 'sk', label="Exp. Upper")
        plt.plot(x_lower[i], C_p_lower[i], 'vk', label="Exp. Lower")
        plt.plot(x_panel, C_p_panel, 'k--', label="Panel Method")
        plt.xlabel("$x/c$")
        plt.ylabel("$C_p$")
        plt.legend()
        plt.gca().invert_yaxis()
        plt.savefig("knight -8 station {0}.pdf".format(i+1))
        plt.savefig("knight -8 station {0}.svg".format(i+1))
        plt.close()