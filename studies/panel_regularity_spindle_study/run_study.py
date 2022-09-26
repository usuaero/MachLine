import sys
import json
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt

sys.path.append("dev/helper_scripts")

from geometry_creator import generate_random_spindle
from paraview.simple import LegacyVTKReader, SaveData

if __name__=="__main__":

    # Parameters
    input_file = "studies/panel_regularity_spindle_study/spindle_input.json"
    mesh_file = "studies/panel_regularity_spindle_study/meshes/random_spindle.vtk"
    result_file = "studies/panel_regularity_spindle_study/results/random_spindle.vtk"
    data_file = "studies/panel_regularity_spindle_study/results/random_spindle_data.csv"
    analytic_data_file = "studies/panel_regularity_spindle_study/method_of_char_from_ehlers.csv"
    M = np.sqrt(2)

    # Generate mesh
    def r_of_x(x):
        return 0.2*x*(1.0-x)
    #generate_random_spindle(mesh_file, 10000, 1.0, r_of_x)

    # Write input
    input_dict = {
        "flow": {
            "freestream_velocity": [1.0, 0.0, 0.0],
            "gamma" : 1.4,
            "freestream_mach_number" : M
        },
        "solver" : {
            #"control_point_offset" : 0.01
            #"matrix_solver" : "GMRES"
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
    spindle_vtk = LegacyVTKReader(FileNames=[result_file])
    SaveData(data_file, proxy=spindle_vtk, FieldAssociation="Cell Data")
    ML_data = np.genfromtxt(data_file, delimiter=',', skip_header=1)
    C_p_ML = ML_data[:,7]
    x_ML = ML_data[:,4]

    # Extract analytic data
    ML_data = np.genfromtxt(analytic_data_file, delimiter=',', skip_header=2)
    C_p_anl = ML_data[:,1]
    x_anl = ML_data[:,0]

    # Plot
    plt.figure()
    plt.plot(x_ML, C_p_ML, 'k.', label='Isentropic')
    plt.plot(x_anl, C_p_anl, 'k', label='Meth. of Char.')
    plt.xlabel("$x$")
    plt.ylabel("$C_P$")
    plt.show()