# This script is to run automate running machline for the Weber and Brebner results

import numpy as np
import json
import subprocess
import time

def mach_iter(AoA_list, Nodes, formulation, freestream):

    if formulation == "source-free":
        formulation_adjusted = "source_free"
    
    else:
        formulation_adjusted = formulation

    filebase = "dev/results/half_wing_swept_45deg/MachLine_Results/"



    for AoA in AoA_list:

        # Modify freestream velocities based on angle of attack
        AoA_deg = float(AoA)*np.pi/180
        x_flow = freestream * np.cos(AoA_deg)
        z_flow = freestream * np.sin(AoA_deg)

        for Node in Nodes:
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
                    "file": "dev/meshes/half_wing_A_meshes/half_wing_A_" + Node + "_nodes.vtk",
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
                    "body_file": filebase +  AoA + "_degrees_AoA/half_wing_A_" + Node + "_nodes_" + AoA + "_deg_AoA.vtk",
                    "wake_file": filebase +  AoA + "_degrees_AoA/half_wing_A_" + Node + "_nodes_" + AoA + "_deg_AoA_wake.vtk",
                    "control_point_file": filebase +  AoA + "_degrees_AoA/half_wing_A_" + Node + "_nodes_" + AoA + "_deg_AoA_" + formulation_adjusted + "_control_points.vtk",
                    "report_file": "dev/report.txt"
                }
            }
            
            # Identify output file location
            filename = AoA + "_deg_angle_of_attack_input.json"
            inputfile = 'dev/half_wing_A_swept_inputs/' + filename
            
            # file_location = "dev/results/half_wing_swept_45deg/test/" + AoA + "_degree_AoA_test_file_" + Node + "_nodes.json"
            with open(inputfile, "w") as output_file:
                json.dump(dict1, output_file, indent=4) 
            
            print("\n***",Node, "node input file saved successfully ***\n")
            
            # Run machline with current input file
            machline_command = "./machline.exe {0}".format(inputfile)
            subprocess.run(machline_command, shell=True)
  

## Main

input_conditions = "dev/results/half_wing_swept_45deg/Swept_half_wing_conditions_input.json"

json_string = open(input_conditions).read()
json_vals = json.loads(json_string)


# Identify values to pass from input conditions file

Nodes_input = json_vals["geometry"]["nodes"]
AoA_list_input = json_vals["geometry"]["AoA list"]
freestream_velocity = json_vals["flow conditions"]["freestream velocity"]

formulation_input = json_vals["solver"]["formulation"]


# Record and print the time required to run MachLine
start_time = time.time()
# Call the machline iterator with the desired inputs
mach_iter(AoA_list_input, Nodes_input, formulation_input, freestream_velocity)    
print("MachLine executed successfully in %s seconds" % "{:.4f}".format(time.time()-start_time))