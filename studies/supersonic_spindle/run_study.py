import json

import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt

from dev.helper_scripts.geometry_creator import generate_random_spindle
from studies.paraview_functions import extract_all_data, get_data_column_from_array
from paraview.simple import LegacyVTKReader, SaveData


RERUN_MACHLINE = True

# Parameters
input_file = "studies/supersonic_spindle/spindle_input.json"
mesh_file = "studies/supersonic_spindle/meshes/random_spindle.vtk"
data_file = "studies/supersonic_spindle/results/random_spindle_data.csv"
analytic_data_file = "studies/supersonic_spindle/method_of_char_from_ehlers.csv"
M = np.sqrt(2)


def run_random_spindle(N):
    """Runs a random spindle with N vertices through MachLine and plots the pressures."""

    # Name results file for this number of vertices
    result_file = "studies/supersonic_spindle/results/random_spindle_{0}.vtk".format(N)

    # Generate mesh
    def r_of_x(x):
        return 0.2*x*(1.0-x)
    generate_random_spindle(mesh_file, N, 1.0, r_of_x)

    # Write input
    input_dict = {
        "flow": {
            "freestream_velocity": [1.0, 0.0, 0.0],
            "gamma" : 1.4,
            "freestream_mach_number" : M
        },
        "solver" : {
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "wake_model" : {
                "wake_present" : False
            }
        },
        "output": {
            "body_file" : result_file,
            "control_point_file" : result_file.replace(".vtk", "_control_points.vtk")
        }
    }
    with open(input_file, 'w') as input_handle:
        json.dump(input_dict, input_handle, indent=4)

    # Run
    sp.run(["./machline.exe", input_file])

    # Extract data from MachLine
    headers, data = extract_all_data(result_file, which_data='cell')
    x_ML = get_data_column_from_array(headers, data, 'centroid:0')
    C_p_ML = get_data_column_from_array(headers, data, 'C_p_ise')

    # Extract analytic data
    ML_data = np.genfromtxt(analytic_data_file, delimiter=',', skip_header=2)
    C_p_anl = ML_data[:,1]
    x_anl = ML_data[:,0]

    # Plot
    plt.figure()
    plt.plot(x_ML, C_p_ML, 'k.', markersize=3, label='Random Mesh')
    plt.plot(x_anl, C_p_anl, 'k', label='Meth. of Char.')
    plt.xlabel("$x$")
    plt.ylabel("$C_P$")
    plt.legend(fontsize=6, title_fontsize=6)
    plt.savefig('studies/supersonic_spindle/plots/M_{0}_{1}.pdf'.format(M, N))
    plt.savefig('studies/supersonic_spindle/plots/M_{0}_{1}.svg'.format(M, N))
    plt.close()


if __name__=="__main__":

    # Numbers of vertices
    Ns = [50, 100, 200, 400, 800]

    for N in Ns:
        run_random_spindle(N)