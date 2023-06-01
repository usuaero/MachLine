import numpy as np
import os
import sys

sys.path.insert(0, './studies')
from case_running_functions import write_input_file, run_quad, cases, quad_labels

RERUN_MACHLINE = True
study_dir = "studies/supersonic_variable_taper_ratio/"

def run_wing_quad(R_T, M, grid, clustered):
    """Runs a case quad for the given Mach number and mesh density"""

    # Storage Locations
    if clustered:
        case_name = "M_{0}_R_T_{1:.0f}_{2}_clustered".format(M,R_T*100,grid)
        mesh_file = study_dir + "meshes/wing_R_T_{0:.0f}_{1}_{2}.vtk".format(R_T*100, grid, "clustered")
    else: 
        case_name = "M_{0}_R_T_{1:.0f}_{2}_non-clustered".format(M,R_T*100,grid)
        mesh_file = study_dir + "meshes/wing_R_T_{0:.0f}_{1}_{2}.vtk".format(R_T*100, grid, "non-clustered")
    if not os.path.exists(study_dir + "results/"):
        os.mkdir(study_dir + "results/")
    results_file = study_dir + "results/"+case_name+".vtk"
    if not os.path.exists(study_dir + "reports/"):
        os.mkdir(study_dir + "reports/")
    report_file = study_dir + "reports/"+case_name+".json"

    input_dict = {
        "flow": {
            "freestream_velocity": [1.0,0.0,0.0],
            "gamma": 1.4,
            "freestream_mach_number": M
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis": "+z",
            "mirror_about": "xy",
            "reference": {
                "area": 1.0
            }
        },
        "solver": {
            "formulation": "morino",
            "sort_system": False
        },
        "post_processing" : {
            "pressure_rules" : {
                "second-order" : True,
                "isentropic" : True,
                "slender-body" : True,
                "linear" : True
            }
        },
        "output" : {
            "body_file" : results_file,
            "report_file" : report_file
        }
    }

    # Dump
    input_file = study_dir + "input.json"
    write_input_file(input_dict, input_file)

    # Run quad
    reports = run_quad(input_file, run=RERUN_MACHLINE)

if __name__ == "__main__":

    # Parameters
    R_Ts = [0,0.25,0.5,0.75,1]
    grids = ["ultra_coarse", "coarse", "medium", "fine"]
    Ms = [1.5, 2, 2.5, 3, 3.5, 4]
    clustering = [True, False]

    # Loop through and run cases
    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):
            for k, clustered in enumerate(clustering):
                for l, R_T in enumerate(R_Ts):
                    run_wing_quad(R_T, M, grid, clustered)

