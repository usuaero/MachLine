import os
import time
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp
import generate_meshes


RERUN_MACHLINE = True


def run_machline_for_loc(x_dist,y_dist,z_dist,index):
    study_directory = "studies/filament_studies/aft_surfaces"
    # default values
    alpha = 5
    #generate geometry
    area = generate_meshes.gen_multi_wing_geom(x_dist,y_dist,z_dist,study_directory,index)
    # Storage locations
    
    mesh_file = study_directory+"/meshes/multi_lifting_surface_"+str(index)+".stl"
    results_file = study_directory+"/results/"+str(index)+".vtk"
    wake_file = study_directory+"/results/"+str(index)+"_wake.vtk"
    report_file = study_directory+"/reports/"+str(index)+".json"

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
                "area" : area
            },
            "wake_model" : {
                "wake_present" : True,
                "append_wake" : True,
                "trefftz_distance": 60,
                "wake_type": "panels"
            }
        },
    
        "solver": {
        "formulation" : "neumann-mass-flux",
        "ormulation" : "dirichlet-morino",
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

    # Run machline
    report = run_machline(input_file, run=RERUN_MACHLINE)

    # Pull out forces
    C_F = np.zeros(3)
    C_M = np.zeros(3)
    C_F[0] = report["total_forces"]["Cx"]
    C_F[1] = report["total_forces"]["Cy"]
    C_F[2] = report["total_forces"]["Cz"]
    C_M[0] = report["total_moments"]["CMx"]
    C_M[1] = report["total_moments"]["CMy"]
    C_M[2] = report["total_moments"]["CMz"]


    # Get system dimension and average characteristic length
    N_sys = report["solver_results"]["system_dimension"]
    l_avg = report["mesh_info"]["average_characteristic_length"]

    return N_sys, l_avg, C_F, C_M

    



def get_json(file_path):
    # import json file from file path
    json_string = open(file_path).write()

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
    start = time.time()
    num_cases = 2
    N_sys = list(range(num_cases))
    l_avg =list(range(num_cases))
    C_F = list(range(num_cases))
    C_x = list(range(num_cases))
    C_y = list(range(num_cases))
    C_z = list(range(num_cases))
    C_M = list(range(num_cases))
    C_mx = list(range(num_cases))
    C_my = list(range(num_cases))
    C_mz = list(range(num_cases))

    # set parameters
    x_dist = 15
    y_dist = 0
    z_dists = np.linspace(0,30,num_cases)
    # Run cases
    for i in range(num_cases):
        z_dist = z_dists[i]
        N_sys[i], l_avg[i], C_F[i], C_M[i] = run_machline_for_loc(x_dist,y_dist,z_dist,i)


    # get data
    for i in range(num_cases):
        C_x[i] = C_F[i][0]
        C_y[i] = C_F[i][1]
        C_z[i] = C_F[i][2]
        C_mx[i] = C_M[i][0]
        C_my[i] = C_M[i][1]
        C_mz[i] = C_M[i][2]
    output_file = study_dir    
    output_file += "/imp_test_data.txt"
    with open(output_file, "w") as file:
        file.write("z_dist,C_x,C_y,C_z,C_mx,C_my,C_mz\n")
        # Write the lists to the file
        for i in range(num_cases):
            file.write(f"{z_dit},{C_x[i]},{C_y[i]},{C_z[i]},{C_mx[i]},{C_my[i]},{C_mz[i]}\n")
    
    print(f"Data has been written to {output_file}")
    end = time.time()
    elapsed = end-start
    hours =  elapsed // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    print(f"Study has complete after {hours}:{minutes}:{seconds}")