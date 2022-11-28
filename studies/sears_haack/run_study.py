from pdb import runcall
import sys
import json
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt

from paraview.simple import LegacyVTKReader, SaveData


def run_case(refinement, M, run_machline=True):

    # Parameters
    input_file = "studies/sears_haack/sears_haack_input.json"
    mesh_file = "studies/sears_haack/meshes/SH_{0}.tri".format(refinement)
    result_file = "studies/sears_haack/results/SH_{0}_M_{1}.vtk".format(refinement, M)
    data_file = "studies/sears_haack/results/SH_{0}_M_{1}_data.csv".format(refinement, M)
    report_file = "studies/sears_haack/reports/SH_{0}_M_{1}.json".format(refinement, M)


    if run_machline:

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
                },
                "reference" : {
                    "area" : 1.675e-3
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
                "body_file" : result_file,
                "report_file" : report_file
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

    # Get report
    with open(report_file, 'r') as report_handle:
        report = json.load(report_handle)

    return x, Cp_ise, Cp_2nd, Cp_lin, Cp_sln, report


if __name__=="__main__":

    # Study parameters
    Ms = [2.0, 3.0]#, 4.0, 6.0, 8.0]
    refinements = ["20_8", "30_11", "40_15", "60_23", "80_30", "120_45", "160_60"]#, "320_120"]

    # Get Stivers data
    data = np.genfromtxt("studies/sears_haack/SH_data_from_stivers.csv", skip_header=1, delimiter=',')

    # Loop through Mach numbers
    N_verts = []
    CDs = []
    for refinement in refinements:
        CDs.append([])
        for i, M in enumerate(Ms):

            try:
                x, Cp_ise, Cp_2nd, Cp_lin, Cp_sln, report = run_case(refinement, M, run_machline=False)
            except:
                continue
            c = np.max(x)
            x_c = x/c
            if i == 0:
                N_verts.append(report["solver_results"]["system_dimension"])
            CDs[-1].append(report["total_forces"]["Cx"])

            # Plot
            N_theta = int(refinement.split('_')[1])
            plt.figure()
            plt.plot(data[:,i*2], data[:,i*2+1], 'k-', label="Meth. of Char.")
            plt.plot(x_c[::N_theta], Cp_ise[::N_theta], 'kv', markersize=3, label="Isentropic")
            plt.plot(x_c[::N_theta], Cp_2nd[::N_theta], 'ks', markersize=3, label="Second-Order")
            plt.plot(x_c[::N_theta], Cp_lin[::N_theta], 'k^', markersize=3, label="Linear")
            plt.plot(x_c[::N_theta], Cp_sln[::N_theta], 'ko', markersize=3, label="Slender-Body")
            plt.xlabel('$x$')
            plt.ylabel('$C_P$')
            if M == 2.0:
                plt.gca().set_ylim((-0.05, 0.4))
            else:
                plt.gca().set_ylim((-0.05, 0.3))
            plt.legend(fontsize=6)
            plt.savefig('studies/sears_haack/plots/SH_{0}_M_{1}.pdf'.format(refinement, M))
            plt.savefig('studies/sears_haack/plots/SH_{0}_M_{1}.svg'.format(refinement, M))
            plt.close()

    # Plot convergence
    CDs = np.array(CDs)
    CD_err = np.abs((CDs[:-1,:] - CDs[-1,:]) / CDs[-1,:])
    plt.figure()
    plt.loglog(N_verts[:-1], CD_err[:,0], color='k', label='2')
    plt.loglog(N_verts[:-1], CD_err[:,1], color='darkgray', label='3')
    plt.xlabel("$N_{vertices}$")
    plt.ylabel("\% Error in $C_D$")
    plt.legend(title="$M_\infty$")
    plt.savefig('studies/sears_haack/plots/SH_convergence.pdf')
    plt.savefig('studies/sears_haack/plots/SH_convergence.svg')
    plt.close()