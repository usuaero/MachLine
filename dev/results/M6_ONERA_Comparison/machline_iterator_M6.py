# This script is to run automate running machline for the Weber and Brebner results

import numpy as np
import json
import subprocess
import time
import multiprocessing as mp
import os
from sys import exit

# Record and print the time required to run MachLine
start_time = time.time()

def mach_iter(AoA, calculation_type, freestream_mach_number):

    # Convert angle of attack to radians
    AoA_rad = float(AoA)*np.pi/180

    # Determine applicible categories for each input file
    speed_of_sound = 340.26 # speed of sound in m/s

    if calculation_type == "pressure_corrections":
        # Adjust freestream settings
        corrected_mach = float(freestream_mach_number)
        isentropic_mach = 0.0
        force_pressure = "prandtl-glauert"

        # Modify freestream velocities based on angle of attack
        x_flow = speed_of_sound * corrected_mach * np.cos(AoA_rad)
        z_flow = speed_of_sound * corrected_mach* np.sin(AoA_rad)

        # Set pressure rules
        pg_rule = True
        kt_rule = True
        lai_rule = True
        incomp_rule = True
        isen_rule = False
        second_rule = False

        
    
    elif calculation_type == "subsonic_calcs":
        # Adjust freestream settings
        corrected_mach = 0.0
        isentropic_mach = float(freestream_mach_number)
        force_pressure = "isentropic"

        # Modify freestream velocities based on angle of attack
        x_flow = speed_of_sound * isentropic_mach * np.cos(AoA_rad)
        z_flow = speed_of_sound * isentropic_mach * np.sin(AoA_rad)
                
        # Set pressure rules
        pg_rule = False
        kt_rule = False
        lai_rule = False
        incomp_rule = False
        isen_rule = True
        second_rule = True

    else:
        print(calculation_type, " is an invalid calculation type. Quitting...")
        exit()




    # Identify filebases used throughout iterator
    filebase = "dev/results/M6_ONERA_Comparison/"
    output_filebase = filebase + "vtk_Results/MachLine/" + calculation_type + f"/MachLine_mach_{freestream_mach_number}_AoA_{AoA}"

    # Rewrite the input files based on angle of attack and mach number
    dict1 = {
        "flow": {
            "freestream_velocity": [x_flow, 0.0, z_flow],
            "gamma": 1.4,
            "freestream_mach_number": isentropic_mach
        },
        "geometry": {
            "file": "dev/meshes/m6_onera.vtk",
            "mirror_about": "xz",
  
            "wake_model": {
                "wake_shedding_angle": 90.0,
                "trefftz_distance": 1000.0,
                "N_panels": 1
            },
            "reference": {
                "area": 1.0
            }
        },
        "solver": {
            "formulation": ["morino"],
            "control_point_offset": 1.1e-05
        },
        "post_processing" : {
            "pressure_for_forces": force_pressure,
            "pressure_rules": {
                "incompressible": incomp_rule,
                "isentropic": isen_rule,
                "second-order": second_rule
        },
        "subsonic_pressure_correction": {
            "correction_mach_number": corrected_mach,
            "prandtl-glauert": pg_rule,
            "karman-tsien": kt_rule,
            "laitone": lai_rule
        }
        },
        "output": {
            "body_file": output_filebase + "_morino.vtk",
            "wake_file": output_filebase + "_morino_wake.vtk",
            "control_point_file": output_filebase + "_morino_control_points.vtk",
            "report_file": "../../report.json"
        }
    }
    
    # Identify output file location
    filename = f"M6_input_file_{AoA}_AoA.json"
    inputfile = filebase + "Input_Files/" + calculation_type + "/" + filename
    
    # file_location = "dev/results/half_wing_swept_45deg/test/" + AoA + "_degree_AoA_test_file_" + Node + "_nodes.json"
    with open(inputfile, "w") as output_file:
        json.dump(dict1, output_file, indent=4) 
    
    # print("\n***",Node, "node input file saved successfully ***\n")
    
    # Run machline with current input file
    subprocess.call(["./machline.exe", inputfile])

## Main


input_conditions = "M6_Onera_input_settings.json"
with open(input_conditions, "r") as json_string:
    json_vals = json.load(json_string)


# Identify values to pass from input conditions file

AoA_list_input = json_vals["geometry"]["AoA list"]
calculations = json_vals["solver"]["calculation_type"]
mach_numbers = json_vals["flow conditions"]["mach_numbers"]

# Identify number of CPU available to work with
# n_processors = mp.cpu_count()
n_processors = 4

Arguments = []

# Change the working directory to the main MachLine directory for execution
os.chdir("../../../")

# Call the machline iterator with the desired inputs
with mp.Pool(n_processors) as pool:
    for calc in calculations:
        for i, AoA in enumerate(AoA_list_input):
            Arguments.append((AoA, calc, mach_numbers[i]))
        pool.starmap(mach_iter, Arguments)


pool.join()

print("MachLine Iterator executed successfully in %s seconds" % "{:.4f}".format(time.time()-start_time))