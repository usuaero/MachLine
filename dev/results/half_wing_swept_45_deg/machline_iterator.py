# This script is to run automate running machline for the Weber and Brebner results

import numpy as np
import json
import subprocess
import time
import multiprocessing as mp
import os

# Record and print the time required to run MachLine
start_time = time.time()

def mach_iter(AoA, Node, formulation, freestream):

    if formulation == "source-free":
        formulation_adjusted = "source_free"
    
    else:
        formulation_adjusted = formulation


    # Modify freestream velocities based on angle of attack
    AoA_rad = float(AoA)*np.pi/180
    x_flow = freestream * np.cos(AoA_rad)
    z_flow = freestream * np.sin(AoA_rad)


    # Identify filebases used throughout iterator
    filebase = "dev/results/half_wing_swept_45_deg/"
    output_filebase = filebase + "MachLine_Results/" +  AoA + "_degrees_AoA/half_wing_A_" + Node + "_nodes_" + AoA + "_deg_AoA_" + formulation_adjusted

    # Rewrite the input files based on angle of attack and node densities
    dict1 = {
        "flow": {
            "freestream_velocity": [
                x_flow,
                0.0,
                z_flow
            ]
        },
        "geometry": {
            "file": filebase + "half_wing_A_meshes/half_wing_A_" + Node + "_nodes.vtk",
            "mirror_about": "xz",
            "singularity_order": {
                "doublet": 1,
                "source": 0
            },
            "wake_model": {
                "wake_shedding_angle": 90.0,
                "trefftz_distance": 10000.0,
                "N_panels": 1
            },
            "reference": {
                "area": 1.0
            }
        },
        "solver": {
            "formulation": formulation,
            "control_point_offset": 1.1e-05
        },
        "post_processing" : {
        },
        "output": {
            "body_file": output_filebase + "_formulation.vtk",
            "wake_file": output_filebase + "_formulation_wake.vtk",
            "control_point_file": output_filebase + "_control_points.vtk",
            "report_file": "../../report.txt"
        }
    }
    
    # Identify output file location
    filename = AoA + "_deg_angle_of_attack_input.json"
    inputfile = filebase + 'half_wing_A_swept_inputs/' + filename
    
    # file_location = "dev/results/half_wing_swept_45deg/test/" + AoA + "_degree_AoA_test_file_" + Node + "_nodes.json"
    with open(inputfile, "w") as output_file:
        json.dump(dict1, output_file, indent=4) 
    
    print("\n***",Node, "node input file saved successfully ***\n")
    
    # Run machline with current input file
    # machline_command = "./machline.exe {0}".format(inputfile)
    subprocess.call(["./machline.exe", inputfile])

## Main


input_conditions = "Swept_half_wing_conditions_input.json"

json_string = open(input_conditions).read()
json_vals = json.loads(json_string)


# Identify values to pass from input conditions file

Nodes_input = json_vals["geometry"]["nodes"]
AoA_list_input = json_vals["geometry"]["AoA list"]
freestream_velocity = json_vals["flow conditions"]["freestream velocity"]

formulation_input = json_vals["solver"]["formulation"]

# Identify number of CPU available to work with
n_processors = mp.cpu_count() * 2

Arguments = []

# Change the working directory to the main MachLine directory for execution
os.chdir("../../../")

# Call the machline iterator with the desired inputs
with mp.Pool(n_processors) as pool:
    for form in formulation_input:
        for AoA in AoA_list_input:
            for node in Nodes_input:
                Arguments.append((AoA, node, form, freestream_velocity))
        pool.starmap(mach_iter, Arguments)


pool.join()

# mach_iter(AoA_list_input, Nodes_input, formulation_input, freestream_velocity)    
print("MachLine Iterator executed successfully in %s seconds" % "{:.4f}".format(time.time()-start_time))