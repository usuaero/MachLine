import json
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt


def run_machline(input_dict, input_filename):
    # Runs MachLine using the provided input dictionary. Will write out to the specified filename

    # Write input
    with open(input_filename, 'w') as input_handle:
        json.dump(input_dict, input_handle)

    # Run MachLine
    sp.run(["./machline.exe", input_filename])

    # Get report
    report_file = input_dict["output"].get("report_file")
    if report_file is not None:
        with open(report_file) as report_handle:
            report = json.load(report_handle)
            return report


if __name__=="__main__":

    # Options
    densities = ["ultra_coarse", "very_coarse", "coarse", "medium"]
    phis = np.radians([30.0, 45.0, 60.0])
    thetas = np.radians([30.0, 45.0, 60.0])
    N = []
    Cx = np.zeros((len(phis), len(thetas), len(densities)))
    Cy = np.zeros((len(phis), len(thetas), len(densities)))
    Cz = np.zeros((len(phis), len(thetas), len(densities)))

    for i, phi in enumerate(phis):
        for j, theta in enumerate(thetas):
            for k, density in enumerate(densities):

                # Initialize input
                result_file = "studies/convergence_study/results/sphere_{0}.vtk".format(density)
                report_file = "studies/convergence_study/reports/sphere_{0}.json".format(density)

                # Lower-order results
                input_dict ={
                    "flow" : {
                        "freestream_velocity" : [np.cos(phi)*np.cos(theta), np.sin(phi)*np.cos(theta), np.sin(theta)]
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

                # Run
                report = run_machline(input_dict, "studies/convergence_study/input.json")

                # Store data
                if i == 0 and j == 0:
                    N.append(report["solver_results"]["system_dimension"])
                Cx[i,j,k] = report["total_forces"]["Cx"]
                Cy[i,j,k] = report["total_forces"]["Cy"]
                Cz[i,j,k] = report["total_forces"]["Cz"]

    # Calculate convergence rates
    orders = []
    for i in range(len(phis)):
        for j in range(len(thetas)):

            # Cx
            err = abs(Cx[i,j,:])
            coefs = np.polyfit(np.log(N), np.log(err), deg=1)
            orders.append(coefs[0])

            # Cy
            err = abs(Cy[i,j,:])
            coefs = np.polyfit(np.log(N), np.log(err), deg=1)
            orders.append(coefs[0])

            # Cz
            err = abs(Cz[i,j,:])
            coefs = np.polyfit(np.log(N), np.log(err), deg=1)
            orders.append(coefs[0])

    avg_order = np.average(orders)
    print("Average order of convergence: ", avg_order)

    # plot
    plt.figure()
    for i in range(len(phis)):
        for j in range(len(thetas)):
            plt.plot(N, np.abs(Cx[i,j,:]), 'k-')
            plt.plot(N, np.abs(Cy[i,j,:]), 'k-')
            plt.plot(N, np.abs(Cz[i,j,:]), 'k-')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$N_{verts}$')
    plt.ylabel('$C_F$')
    plt.show()