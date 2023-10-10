import os
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp
import generate_meshes


RERUN_MACHLINE = True


def run_machline_for_loc(x_dist,y_dist_z_dist):
    # default values
    alpha = 5
    # Storage locations
    case_name = "{0}_{1}_{2}".format(x_dist,y_dist_z_dist)
    mesh_file = "studies/filament_studies/aft_surfaces/meshes/multi_lifting_surface_"+case_name+".stl"
    results_file = "studies/filament_studies/aft_surfaces/results/"+case_name+".vtk"
    wake_file = "studies/filament_studies/aft_surfaces/results/"+case_name+"_wake.vtk"
    report_file = "studies/filament_studies/aft_surfaces/reports/"+case_name+".json"

    # create input file
    input_dict = {
        "flow": {
            "freestream_velocity": [np.cos(np.radians(alpha)), 0.0, np.sin(np.radians(alpha))],
            "freestream_mach_number" : 0.25
            
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "reference" : {
                "area" : 50.0
            },
            "wake_model" : {
                "wake_present" : True,
                "append_wake" : True,
                "trefftz_distance": trefftz,
                "wake_type": "panel"
            }
        },
    
        "solver": {
        "ormulation" : "neumann-mass-flux",
        "formulation" : "dirichlet-morino",
        "matrix_solver" : "GMRES",
        "write_A_and_b" : False
        },
        "post_processing" : {
            "pressure_rules" : {
                "ncompressible" : True
            }
        },
        "output" : {
            "body_file" : results_file,
            "wake_file" : wake_file,
            "report_file" : report_file
        }
    }
    # Dump
    input_file = "studies/filament_studies/aft_surfaces/input.json"
    write_input_file(input_dict, input_file)

    # Run quad
    report = run_machline(input_file, run=RERUN_MACHLINE)

    # Pull out forces
    C_F = np.zeros(3)
    C_F[0] = report["total_forces"]["Cx"]
    C_F[1] = report["total_forces"]["Cy"]
    C_F[2] = report["total_forces"]["Cz"]

    # Get system dimension and average characteristic length
    N_sys = report["solver_results"]["system_dimension"]
    l_avg = report["mesh_info"]["average_characteristic_length"]

    return N_sys, l_avg, C_F

    



def get_json(file_path):
    # import json file from file path
    json_string = open(file_path).read()

    # save to vals dictionary
    input_dict = json.loads(json_string)
    
    return input_dict

def write_input_file(input_dict, input_filename):
    """Writes the given input dict to the given file location."""

    with open(input_filename, 'w') as input_handle:
        json.dump(input_dict, input_handle, indent=4)

def run_machline(input_filename, delete_input=True, run=True):
    """Runs MachLine with the given input and returns the report if MachLine generated one."""

    # Run
    if run:
        sp.run(["./machline.exe", input_filename])

    # Get report
    with open(input_filename, 'r') as input_handle:
        input_dict = json.load(input_handle)
    report_file = input_dict["output"].get("report_file")
    if report_file is not None:
        try:
            with open(report_file) as report_handle:
                report = json.load(report_handle)
        except:
            report = None
    else:
        report = None

    # Delete input
    if delete_input:
        os.remove(input_filename)

    return report


if __name__=="__main__":
    # declare varaibles and inputs
    num_cases = 30
    N_sys = list(range(num_cases))
    l_avg =list(range(num_cases))
    C_F = list(range(num_cases))
    C_x = list(range(num_cases))
    C_y = list(range(num_cases))
    C_z = list(range(num_cases))

    # Run cases
    trefftz_distances = np.linspace(1,500,num_cases)
    for i in range(num_cases):
        N_sys[i], l_avg[i], C_F[i] = run_machline_for_z_dist(trefftz_distances[i])


    # get data
    for i in range(num_cases):
        C_x[i] = C_F[i][0]
        C_y[i] = C_F[i][1]
        C_z[i] = C_F[i][2]
        
    #  plot
    plt.figure()
    plt.plot(trefftz_distances,C_x, label="C_x")
    plt.plot(trefftz_distances,C_y,label="C_y")
    plt.plot(trefftz_distances,C_z,label="C_z")
    
    plt.xlabel("Trefftz Distance")
    plt.ylabel("$C_F$")
    plt.legend()
    plt.show()

