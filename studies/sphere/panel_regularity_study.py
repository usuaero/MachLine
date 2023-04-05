import json

import numpy as np
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt

from studies.paraview_functions import extract_all_data, get_data_column_from_array
from studies.case_running_functions import write_input_file, run_quad, cases, line_styles


RERUN_MACHLINE = True


def run_sphere_comparison(grid, sample, run_machline=True):
    # Runs the comparison of the random sphere

    # Storage locations
    case_name = "random_sphere_{0}_sample_{1}".format(grid, sample)
    plot_dir = "studies/sphere/plots/"
    mesh_file = "studies/sphere/meshes/"+case_name+".vtk"
    body_file = "studies/sphere/results/"+case_name+".vtk"
    report_file = "studies/sphere/reports/"+case_name+".json"

    # Declare MachLine input
    input_dict = {
        "flow" : {
            "freestream_velocity": [10.0, 0.0, 0.0]
        },
        "geometry" : {
            "file" : mesh_file,
            "spanwise_axis" : "+y",
            "wake_model" : {
                "wake_present" : False
            },
            "reference" : {
                "area" : np.pi
            }
        },
        "solver": {
            "formulation": "morino",
            "control_point_offset": 1.1e-6,
            "matrix_solver" : "BJAC",
            "relaxation" : 0.9,
            "tolerance" : 1.1e-11,
            "block_size" : 400
        },
        "post_processing" : {
        },
        "output" : {
            "body_file" : body_file,
            "report_file" : report_file
        }
    }

    # Dump
    input_file = "studies/sphere/sphere_input.json"
    write_input_file(input_dict, input_file)

    # Run
    reports = run_quad(input_file, run=RERUN_MACHLINE)
    
    # Extract lift, drag, and pressures
    Cx = []
    Cy = []
    Cz = []
    C_P = []
    x = []
    for i, report in enumerate(reports):

        # Force coefficients
        Cx.append(report["total_forces"]["Cx"])
        Cy.append(report["total_forces"]["Cy"])
        Cz.append(report["total_forces"]["Cz"])

        # Get pressures
        column_headers, data = extract_all_data(report["input"]["output"]["body_file"], which_data='cell')
        C_P.append(get_data_column_from_array(column_headers, data, "C_p_inc"))
        x.append(get_data_column_from_array(column_headers, data, "centroid:0"))

        # Plot data from MachLine
        plt.figure()
        theta = np.degrees(np.arccos(x[i]))
        theta_a = np.linspace(0.0, np.pi, 100)
        plt.plot(theta, C_P[i], 'k.', markersize=1, label=cases[i])
        plt.plot(np.degrees(theta_a), 1.0-2.25*np.sin(theta_a)**2, 'k-', label='Analytic')
        plt.xlabel('$\\theta [^\circ]$')
        plt.ylabel('$C_p$')
        plt.ylim(bottom=-2.5)
        plt.legend(fontsize=6, title_fontsize=6)
        plt.savefig(plot_dir+case_name+"_{0}.pdf".format(cases[i]))
        plt.savefig(plot_dir+case_name+"_{0}.svg".format(cases[i]))
        plt.close()

    return Cx, Cy, Cz


if __name__=="__main__":

    Ns = [125, 250, 500, 1000, 2000]
    grids = ["ultra_coarse", "coarse", "medium", "fine", "ultra_fine"]
    samples = [x for x in range(10)]

    C_f = np.zeros((len(grids), len(samples), 4, 3))

    # Plot each sample
    for i, grid in enumerate(grids):
        for j, sample in enumerate(samples):

                C_f[i,j,:,0], C_f[i,j,:,1], C_f[i,j,:,2] = run_sphere_comparison(grid, sample, run_machline=False)

    ## Plot convergence
    #avg = np.average(C_f, axis=1)
    #std_dev = np.std(C_f, axis=1)
    #plot_dir = "studies/sphere/plots/"
    #plt.figure()
    #plt.plot(Ns, np.abs(avg[:,0]), 'ks-', label="$C_x$")
    #plt.plot(Ns, np.abs(avg[:,1]), 'ko--', label="$C_y$")
    #plt.plot(Ns, np.abs(avg[:,2]), 'kv-.', label="$C_z$")
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.xlabel("$N_{verts}$")
    #plt.ylabel("Force Coefficient")
    #plt.legend()
    #plt.savefig("studies/sphere/plots/convergence.pdf")
    #plt.savefig("studies/sphere/plots/convergence.svg")