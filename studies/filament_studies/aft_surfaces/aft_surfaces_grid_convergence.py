import os
import time
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp
import generate_meshes


RERUN_MACHLINE = True


def run_machline_for_comp(wake_type,freestream_mach,append_wake,formulation,study_directory,both_wing):

    
    # default values
    alpha = 5
    #generate geometry
    # area = generate_meshes.gen_multi_wing_geom(x_dist,y_dist,z_dist,study_directory,index)
    # Storage locations
    area = 1

    caseName = wake_type+"_"+str(freestream_mach).replace(".","_")+"_"+str(append_wake)+"_"+formulation+"_"
    
    if both_wing:
        mesh_file = study_directory+"/meshes/bothWings.stl"
        caseName += "bothWings"
    else:
        mesh_file = study_directory+"/meshes/oneWing.stl"
        caseName += "oneWing"
    
    results_file = study_directory+"/results/"+caseName+".vtk"
    wake_file = study_directory+"/results/"+caseName+"wake.vtk"
    control_pts_file = study_directory+"/results/"+caseName+"_control_pts.vtk"
    report_file = study_directory+"/reports/"+caseName+".json"

    # create input file
    input_dict = {
        "flow": {
            "freestream_velocity": [np.cos(np.radians(alpha)), 0.0, np.sin(np.radians(alpha))],
            "freestream_mach_number" : freestream_mach
            
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "reference" : {
                "area" : area
            },
            "wake_model" : {
                "wake_present" : True,
                "append_wake" : append_wake,
                "trefftz_distance": 60,
                "wake_type": wake_type
            }
        },
    
        "solver": {
        "formulation" : formulation,
        "xmatrix_solver" : "GMRES",
        "matrix_solver" : "HHLS",
        "write_A_and_b" : False,
        "sort_system" : True
        },
        "post_processing" : {
            "pressure_rules" : {
                "ncompressible" : True
            }
        },
        "output" : {
            "body_file" : results_file,
            "wake_file" : wake_file,
            "control_point_file" : control_pts_file,
            "report_file" : report_file,
            "verbose": True
        }
    }
    # Dump
    input_file = "studies/filament_studies/aft_surfaces/input.json"
    write_input_file(input_dict, input_file)

    # Run machline
    report = run_machline(input_file, run=RERUN_MACHLINE)

    # Pull out forces
    c_F = np.zeros(3)
    c_M = np.zeros(3)
    c_F[0] = report["total_forces"]["Cx"]
    c_F[1] = report["total_forces"]["Cy"]
    c_F[2] = report["total_forces"]["Cz"]
    c_M[0] = report["total_moments"]["CMx"]
    c_M[1] = report["total_moments"]["CMy"]
    c_M[2] = report["total_moments"]["CMz"]


    # Get system dimension and average characteristic length
    N_sys = report["solver_results"]["system_dimension"]
    l_avg = report["mesh_info"]["average_characteristic_length"]

    return N_sys, l_avg, c_F, c_M

    



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
    study_directory = "studies/filament_studies/aft_surfaces"
    start = time.time()
    

    # set parameters
    wake_types = ["panels","filaments","dirichlet","noneD","noneN"]
    formulations = ["neumann-mass-flux","dirichlet-morino"]
    freestream_machs = [2.0]
    num_cases = len(freestream_machs)*len(wake_types)

    N_sys = np.zeros(num_cases)
    l_avg = np.zeros(num_cases) 
    c_F =  np.zeros((num_cases,3))
    c_M =  np.zeros((num_cases,3))
    c_FswList = np.zeros((num_cases,3))
    c_MswList = np.zeros((num_cases,3))
    c_FbwList = np.zeros((num_cases,3))
    c_MbwList = np.zeros((num_cases,3))
    index = 0

    print("Running aft surfaces comparison study")
    # Run cases

    
        
# Run cases

    for i in range(len(freestream_machs)):
        freestream_mach = freestream_machs[i]
        for j in range(len(wake_types)):
            # set parameters
            if wake_types[j] == "noneD" or wake_types[j] == "noneN":
                append_wake = False
                wake_type = "panels"
                if wake_types[j] == "noneD":
                    formulation = formulations[1]
                else:
                    formulation = formulations[0]
            else:
                append_wake = True
                if wake_types[j] == "panels":
                    formulation = formulations[0]
                    wake_type = wake_types[j]
                elif wake_types[j] == "filaments":
                    formulation = formulations[0]
                    wake_type = wake_types[j]
                else:
                    formulation = formulations[1]
                    wake_type = "panels"
            
            
            N_sys, l_avg, c_F_sw, c_M_sw = run_machline_for_comp(wake_type,freestream_mach,append_wake,formulation,study_directory,False)
            N_sys, l_avg, c_F_bw, c_M_bw = run_machline_for_comp(wake_type,freestream_mach,append_wake,formulation,study_directory,True)
            c_F[index] = c_F_bw - c_F_sw
            c_M[index] = c_M_bw - c_M_sw
            c_FswList[index] = c_F_sw
            c_MswList[index] = c_M_sw
            c_FbwList[index] = c_F_bw
            c_MbwList[index] = c_M_bw
            index += 1
    
    index = 0
    for i in range(len(freestream_machs)):
        for j in range(len(wake_types)):
            wake_type = wake_types[j]
            freestream_mach = freestream_machs[i]
            print(f"Results for {wake_type} wake with freestream mach {freestream_mach}")              
            print("c_F = ",c_F[index])
            print("c_M = ",c_M[index])
            print("c_Fsw = ",c_FswList[index])
            print("c_Msw = ",c_MswList[index])
            print("c_Fbw = ",c_FbwList[index])
            print("c_Mbw = ",c_MbwList[index])
            index += 1
    
    with open("studies/filament_studies/aft_surfaces/results.txt",'w') as results_handle:
        results_handle.write("Results for aft surfaces comparison study\n")
        index = 0
        for i in range(len(freestream_machs)):
            for j in range(len(wake_types)):
                result = "Results for "+wake_types[j]+" wake with freestream mach "+str(freestream_machs[i])+"\n"
                result += "c_F,"+str(c_F[index])+",c_M,"+str(c_M[index])+",c_Fsw,"+str(c_FswList[index])+",c_Msw,"+str(c_MswList[index])+",c_Fbw,"+str(c_FbwList[index])+",c_Mbw,"+str(c_MbwList[index])+"\n"
                results_handle.write(result)
                index += 1

    end = time.time()
    seconds = end-start
    hours =  seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    print(f"Study has complete after {int(hours)}:{int(minutes)}:{int(seconds)}")


