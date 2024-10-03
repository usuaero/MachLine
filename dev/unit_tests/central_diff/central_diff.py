import os
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp
import time
import sys



RERUN_MACHLINE = True


def run_machline_central_diff(study_directory, cp_offset, perturb_point, point_index, xyz_index, step, formulation):
    
    # print("Running MachLine for cp_offset = {0}".format(cp_offset))
    # Storage locations
    cp_offset_string = str(cp_offset).replace(".","_") + "_" + formulation.replace("-","_")
    case_name = "one_run_in_central_diff".format(cp_offset_string)
    mesh_file = study_directory+"/meshes/test_mesh_11.stl"
    results_file = study_directory+"/results/"+case_name+".vtk"
    # wake_file = study_directory+"/results/"+case_name+"_wake.vtk"
    report_file = study_directory+"/reports/"+case_name+".json"
    control_point_file = study_directory+"/results/"+case_name+"_control_points.vtk"
    

    # create input file
    input_dict = {
        "flow": {
            "freestream_velocity": [1.0, 0.0, 0.1],
            "freestream_mach_number" : 2.0
        },
        "geometry": {
            "file": mesh_file,
            "wake_model" : {
                "wake_present" : True, 
                "wake_appended" : True
            },
            "perturbation" : {
                "perturb_point" : perturb_point,
                "point_index" : point_index,
                "xyz_index" : xyz_index,
                "step" : step
            },
        },
    
        "solver": {
            "formulation" : formulation,
            "control_point_offset": cp_offset,
            "matrix_solver" : "GMRES"
        },
        "post_processing" : {
            "pressure_rules" : {
                "isentropic" : True
            }
        },
        "output" : {
            "verbose" : False,
            "report_file" : report_file
        }
    }
    # Dump
    input_file = study_directory+"/input.json"
    write_input_file(input_dict, input_file)

    # Run MachLine
    report = run_machline(input_file, run=RERUN_MACHLINE)

    # Pull forces
    CF = np.zeros(3)
    CF[0] = report["total_forces"]["Cx"]
    CF[1] = report["total_forces"]["Cy"]
    CF[2] = report["total_forces"]["Cz"]

    # delete report
    # os.remove(report)  # doesnt work

    return CF

    



def get_json(file_path):
    # import json file from file path
    json_string = open(file_path).read()

    # save to vals dictionary
    input_dict = json.loads(json_string)
    
    return input_dict

def write_input_file(input_dict, input_filename):
    """Writes the given input dict to the given file location."""

    with open(input_filename, 'w') as input_handle:
        json.dump(input_dict, input_handle, indent=4,)

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
    file_path = "dev/unit_tests/central_diff/reports/one_run_in_central_diff.json"
    os.remove(file_path)

    return report


if __name__=="__main__":

    # declare varaibles and inputs
    tStart = time.time()
    alpha = np.arctan2(0.1,1.0)
    # num_cases = 1000
    # mach = 2.0
    # wake_present = True
    # wake_appended = True
    formulation = "dirichlet-source-free"
    calc_adjoint = False
    cp_offset = 1.0e-4

    # number of mesh points/vertices
    points = 1190
    step = 1.0e-6

    study_directory = "dev/unit_tests/central_diff"
    
    d_CFx = np.zeros(points*3)
    d_CFy = np.zeros(points*3) 
    d_CFz = np.zeros(points*3)
    
    

    # init plots for values of d_CFx, d_CFy, d_CFz
    figx, ax1 = plt.subplots(figsize=(10,6))
    figy, ax2 = plt.subplots(figsize=(10,6))
    figz, ax3 = plt.subplots(figsize=(10,6))

    perturb_point = True
    point_index = 1
    xyz_index = 1
    run_count = 0

    # xyz loop
    for j in range(1,4):
        
        # vertex loop
        for k in range(1, points + 1):
            
            i = (k + (j-1)*points) -1
            point_index = k
            xyz_index = j

            # print("\nvertex ", k, "    xyz = ",j,"      python index = ",i)
            
            # perturb up
            CF_up = run_machline_central_diff(study_directory, cp_offset, perturb_point, point_index, xyz_index, step, formulation)
            run_count += 1

            # perturb down
            CF_down = run_machline_central_diff(study_directory, cp_offset, perturb_point, point_index, xyz_index, -step, formulation)
            run_count += 1


            # calculate sensitivity wrt this design variable and store
            d_CFx[i] = (CF_up[0] - CF_down[0]) / (2.*step)  
            d_CFy[i] = (CF_up[1] - CF_down[1]) / (2.*step)  
            d_CFz[i] = (CF_up[2] - CF_down[2]) / (2.*step)  

            

        # calc norms
        d_CFx_norm = np.sqrt(np.sum(d_CFx[:] * d_CFx[:]))
        d_CFy_norm = np.sqrt(np.sum(d_CFy[:] * d_CFy[:]))
        d_CFz_norm = np.sqrt(np.sum(d_CFz[:] * d_CFz[:]))

        print("Norm of d_CFx = ", d_CFx_norm)
        print("Norm of d_CFy = ", d_CFy_norm)
        print("Norm of d_CFz = ", d_CFz_norm)
    


    tEnd = time.time()
    print("Elapsed time is {0} seconds".format(tEnd-tStart))

    print("Machline ran " + str(run_count) + " times")
    
    plt.show()
    
    
    