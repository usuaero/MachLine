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
    N = []
    Cx_lower = []
    Cy_lower = []
    Cz_lower = []
    Cx_higher = []
    Cy_higher = []
    Cz_higher = []

    for i, density in enumerate(densities):

        # Initialize input
        result_file = "studies/convergence_study/results/sphere_{0}.vtk".format(density)
        report_file = "studies/convergence_study/reports/sphere_{0}.json".format(density)

        # Lower-order results
        input_dict ={
            "flow" : {
                "freestream_velocity" : [1.0, 1.0, 1.3]
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
        N.append(report["solver_results"]["system_dimension"])
        Cx_lower.append(report["total_forces"]["Cx"])
        Cy_lower.append(report["total_forces"]["Cy"])
        Cz_lower.append(report["total_forces"]["Cz"])

        # Higher-order results
        input_dict ={
            "flow" : {
                "freestream_velocity" : [1.0, 1.0, 1.3]
            },
            "geometry" : {
                "file" : "studies/convergence_study/meshes/sphere_{0}.stl".format(density),
                "singularity_order" : "higher",
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
        Cx_higher.append(report["total_forces"]["Cx"])
        Cy_higher.append(report["total_forces"]["Cy"])
        Cz_higher.append(report["total_forces"]["Cz"])

    # plot
    plt.figure()
    plt.plot(N, np.abs(Cx_lower), 'k-', label='Lower-Order')
    plt.plot(N, np.abs(Cy_lower), 'k-')
    plt.plot(N, np.abs(Cz_lower), 'k-')
    plt.plot(N, np.abs(Cx_higher), 'k--', label='Higher-Order')
    plt.plot(N, np.abs(Cy_higher), 'k--')
    plt.plot(N, np.abs(Cz_higher), 'k--')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$N_{verts}$')
    plt.ylabel('$C_F$')
    plt.legend()
    plt.show()