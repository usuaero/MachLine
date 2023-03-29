import json
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt
from dev.helper_scripts.run_case_quad import *


if __name__=="__main__":

    # Options
    input_filename = "studies/convergence_study/input.json"
    densities = ["ultra_coarse", "very_coarse", "coarse", "medium"]
    psis = np.radians([30.0, 45.0, 60.0])
    thetas = np.radians([30.0, 45.0, 60.0])
    N = []
    l_avg = []
    Cx = np.zeros((len(psis), len(thetas), len(densities)))
    Cy = np.zeros((len(psis), len(thetas), len(densities)))
    Cz = np.zeros((len(psis), len(thetas), len(densities)))

    for i, psi in enumerate(psis):
        for j, theta in enumerate(thetas):
            for k, density in enumerate(densities):

                # Initialize input
                result_file = "studies/convergence_study/results/sphere_{0}.vtk".format(density)
                report_file = "studies/convergence_study/reports/sphere_{0}.json".format(density)

                # Lower-order results
                input_dict ={
                    "flow" : {
                        "freestream_velocity" : [np.cos(psi)*np.cos(theta), np.sin(psi)*np.cos(theta), np.sin(theta)]
                    },
                    "geometry" : {
                        "file" : "studies/convergence_study/meshes/sphere_{0}.stl".format(density),
                        "spanwise_axis" : "+y"
                    },
                    "solver" : {
                    },
                    "post_processing" : {
                    },
                    "output" : {
                        "body_file" : result_file,
                        "report_file" : report_file
                    }
                }

                # Write to file
                with open(input_filename, 'w') as input_handle:
                    json.dump(input_dict, input_handle, indent=4)

                # Run
                report = run_machline(input_filename)

                # Store data
                if i == 0 and j == 0:
                    N.append(report["solver_results"]["system_dimension"])
                    l_avg.append(report["mesh_info"]["average_characteristic_length"])
                Cx[i,j,k] = report["total_forces"]["Cx"]
                Cy[i,j,k] = report["total_forces"]["Cy"]
                Cz[i,j,k] = report["total_forces"]["Cz"]

    # Calculate convergence rates
    orders = []
    for i in range(len(psis)):
        for j in range(len(thetas)):

            # Cx
            err = abs(Cx[i,j,:])
            coefs = np.polyfit(np.log(l_avg), np.log(err), deg=1)
            orders.append(coefs[0])

            # Cy
            err = abs(Cy[i,j,:])
            coefs = np.polyfit(np.log(l_avg), np.log(err), deg=1)
            orders.append(coefs[0])

            # Cz
            err = abs(Cz[i,j,:])
            coefs = np.polyfit(np.log(l_avg), np.log(err), deg=1)
            orders.append(coefs[0])

    # Report average order
    avg_order = np.average(orders)
    print()
    print("Average order of convergence: ", avg_order)
    print("Std dev: ", np.std(orders))

    # plot
    plt.figure()
    for i in range(len(psis)):
        for j in range(len(thetas)):
            plt.plot(l_avg, np.abs(Cx[i,j,:]), 'k-')
            plt.plot(l_avg, np.abs(Cy[i,j,:]), 'k-')
            plt.plot(l_avg, np.abs(Cz[i,j,:]), 'k-')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$l_{avg}$')
    plt.ylabel('$C_F$')
    plt.show()