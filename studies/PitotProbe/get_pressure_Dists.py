import os
import time
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp


RERUN_MACHLINE = True


def run_machline_for_ab(a,b,study_directory):
    x = -np.cos(a)*np.cos(b)/(1-np.sin(a)**2*np.sin(b)**2)**0.5
    y = -np.cos(a)*np.sin(b)/(1-np.sin(a)**2*np.sin(b)**2)**0.5
    z = -np.sin(a)*np.cos(b)/(1-np.sin(a)**2*np.sin(b)**2)**0.5
    
    mesh_file = study_directory+"/meshes/EURP_PitotProbeOpenVSPc.stl"
    results_file = study_directory+"/results/EURP_PitotProbeOpenVSP.vtk"
    report_file = study_directory+"/reports/EURP_PitotProbeOpenVSP.json"

    # create input file
    input_dict = {
        "flow": {
            "freestream_velocity": [x,y,z],
            "freestream_mach_number" : 0
            
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "reference" : {
                "area" : np.pi
            },
            "wake_model" : {
                "wake_present" : False,
                "append_wake" : True,
                "trefftz_distance": 60
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
            "report_file" : report_file,
            "offbody_points": {
                "points_file": study_directory+"/port_locations.csv",
                "output_file": study_directory+"/port_pressures.csv"
            }
        }
    }
    # Dump
    input_file = study_directory + "/input.json"
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
    study_directory = "studies/PitotProbe"
    start = time.time()
    num_cases = 2
    N_sys = list(range(num_cases))
    l_avg =list(range(num_cases))
    C_F = list(range(num_cases))
    C_M = list(range(num_cases))
    
    # set parameters
    a = np.linspace(0,25,num_cases)
    b = np.linspace(0,25,num_cases)
    
    # Run cases
    index = 0
    for i in range(num_cases):
        for j in range(num_cases):
            N_sys[index], l_avg[index], C_F[index], C_M[index] = run_machline_for_ab(a[i],b[j],study_directory)
            index +=1
                
    end = time.time()
    seconds = end-start
    hours =  seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    print(f"Study has complete after {int(hours)}:{int(minutes)}:{int(seconds)}")
    