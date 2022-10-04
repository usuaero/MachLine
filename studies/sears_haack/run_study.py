from pdb import runcall
import sys
import json
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt

from paraview.simple import LegacyVTKReader, SaveData


def run_case(M):

    # Parameters
    input_file = "studies/sears_haack/sears_haack_input.json"
    mesh_file = "studies/sears_haack/meshes/SH_80_30.tri"
    result_file = "studies/sears_haack/results/SH_80_30_M_{0}.vtk".format(M)
    data_file = "studies/sears_haack/results/SH_80_30_M_{0}_data.csv".format(M)

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
        "post_processing" : {
            "pressure_rules" : {
                "isentropic" : True,
                "second-order" : True,
                "slender-body" : True,
                "linear" : True
            }
        },
        "output": {
            "body_file" : result_file
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
    x = ML_data[:,4]
    Cp_ise = ML_data[:,7]
    Cp_2nd = ML_data[:,8]
    Cp_lin = ML_data[:,9]
    Cp_sln = ML_data[:,10]

    return x, Cp_ise, Cp_2nd, Cp_lin, Cp_sln

if __name__=="__main__":

    # Study parameters
    Ms = [2.0, 3.0]#, 4.0, 6.0, 8.0]

    # Get Stivers data
    data = np.genfromtxt("studies/sears_haack/SH_data_from_stivers.csv", skip_header=1, delimiter=',')

    # Loop through Mach numbers
    xs = []
    Cps = []
    for i, M in enumerate(Ms):

        try:
            x, Cp_ise, Cp_2nd, Cp_lin, Cp_sln = run_case(M)
        except:
            continue
        c = np.max(x)
        x_c = x/c

        # Plot
        plt.figure()
        plt.plot(data[:,i*2], data[:,i*2+1], 'k-', label="Meth. of Char.")
        plt.plot(x_c[::30], Cp_ise[::30], 'kv', markersize=3, label="Isentropic")
        plt.plot(x_c[::30], Cp_2nd[::30], 'ks', markersize=3, label="Second-Order")
        plt.plot(x_c[::30], Cp_lin[::30], 'k^', markersize=3, label="Linear")
        plt.plot(x_c[::30], Cp_sln[::30], 'ko', markersize=3, label="Slender-Body")
        plt.xlabel('$x$')
        plt.ylabel('$C_P$')
        if M == 2.0:
            plt.gca().set_ylim((-0.05, 0.4))
        else:
            plt.gca().set_ylim((-0.05, 0.3))
        plt.legend(fontsize=6)
        plt.show()