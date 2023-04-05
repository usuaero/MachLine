import json

import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, write_input_file, cases, line_styles


RERUN_MACHLINE = True


if __name__=="__main__":

    # Sphere
    input_dict = {
        "flow" : {
            "freestream_velocity" : [1.0, 1.0, 1.0]
        },
        "geometry" : {
            "file" : "studies/sphere/meshes/sphere_medium.stl",
            "wake_model" : {
                "wake_present" : False
            },
            "reference" : {
                "area" : np.pi
            }
        },
        "solver" : {
        },
        "post_processing" : {
        },
        "output" : {
        }
    }

    # Set up offsets
    N_offsets = 10
    offsets = np.logspace(-11.999, 0.001, N_offsets)
    Cx = np.zeros((4,N_offsets))
    Cy = np.zeros((4,N_offsets))
    Cz = np.zeros((4,N_offsets))

    input_file = "studies/sphere/input.json"

    for j, offset in enumerate(offsets):

        # Set offset
        input_dict["solver"]["control_point_offset"] = offset
        input_dict["output"]["report_file"] = "studies/sphere/reports/cp_offset_{0}.json".format(round(np.log10(offset)))

        # Write
        write_input_file(input_dict, input_file)

        # Run
        reports = run_quad(input_file, run=RERUN_MACHLINE)
        
        # Load results
        for i, report in enumerate(reports):
            try:
                Cx[i,j] = report["total_forces"]["Cx"]
                Cy[i,j] = report["total_forces"]["Cy"]
                Cz[i,j] = report["total_forces"]["Cz"]
            except:
                Cx[i,:] = np.nan
                Cy[i,:] = np.nan
                Cz[i,:] = np.nan

    # Plot
    plt.figure()
    for i, case in enumerate(cases):
        plt.plot(offsets, np.abs(Cx[i,:]), line_styles[i], label=case)
        plt.plot(offsets, np.abs(Cy[i,:]), line_styles[i])
        plt.plot(offsets, np.abs(Cz[i,:]), line_styles[i])
    plt.xlabel('$k_1$')
    plt.ylabel('$C_F$')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()
    print(Cy)