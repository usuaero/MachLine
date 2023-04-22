import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, cases, line_styles, write_input_file


RERUN_MACHLINE = True
study_dir = "studies/nacelles/"
plot_dir = study_dir + "plots/supersonic/"


def run_study_at_mach(M):
    # Runs the study for the given Mach number

    # options
    grids = ["coarse", "medium", "fine", "ultra_fine"]

    # Loop through grid densities
    for i, grid in enumerate(grids):

        # Files
        mesh_file = study_dir + "meshes/supersonic_{0}.stl".format(grid)
        result_file = study_dir + "results/supersonic_M_{0}_{1}.vtk".format(M, grid)
        control_point_file = study_dir + "results/supersonic_M_{0}_{1}_control_points.vtk".format(M, grid)
        input_file = study_dir + "input.json"
        report_file = study_dir + "reports/supersonic_M_{0}_{1}.json".format(M, grid)

        # Declare input
        input_dict = {
            "flow": {
                "freestream_velocity": [1.0, 0.0, 0.0],
                "gamma" : 1.4,
                "freestream_mach_number" : M
            },
            "geometry": {
                "file": mesh_file,
                "spanwise_axis" : "+y",
                "max_continuity_angle" : 5.0,
                "wake_model": {
                    "append_wake" : False,
                },
                "reference": {
                }
            },
            "solver": {
            },
            "post_processing" : {
                "pressure_rules" : {
                    "isentropic" : True
                }
            },
            "output" : {
                "body_file" : result_file,
                "control_point_file" : control_point_file,
                "report_file" : report_file
            }
        }

        # Write input
        write_input_file(input_dict, input_file)

        # Run
        reports = run_quad(input_file, run=RERUN_MACHLINE)

if __name__=="__main__":

    # Mach numbers
    Ms = [1.5]#, 2.0, 2.5]

    for M in Ms:
        run_study_at_mach(M)