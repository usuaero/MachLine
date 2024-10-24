import os
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp
import time
import sys
import csv
import pandas as pd
import pickle
from write_vtk import read_vtk_file, write_vtk_file, add_vector_data_to_vtk
import concurrent.futures 
from threading import Lock, local


RERUN_MACHLINE = True


def run_machline_for_cp_offset(cp_offset,study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, mesh_name):
    
    print("Running MachLine for cp_offset = {0}".format(cp_offset))
    # Storage locations
    cp_offset_string = str(cp_offset).replace(".","_") + "_" + formulation.replace("-","_")
    case_name = "cp_offset_{0}".format(cp_offset_string)
    mesh_file = study_directory+"/meshes/"+mesh_name+".stl"

    if (calc_adjoint):
        report_file = study_directory + "/reports/adjoint/"+ f'{cp_offset:.2e}' + ".json"
        body_file = study_directory + "/vtk_files/adjoint_"+mesh_name+"_cp_"+f'{cp_offset:.2e}'+".vtk"
    else:
        report_file = study_directory+"/reports/" + str(point_index)  + "_" + str(xyz_index) + "_" + f'{step: .2e}' + f'{cp_offset:.2e}' + ".json"
        body_file = "none"
    

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
            "adjoint_sensitivities" : {
                "calc_adjoint" : calc_adjoint
            },
            "perturbation" : {
                "perturb_point" : perturb_point,
                "point_index" : point_index,
                "xyz_index" : xyz_index,
                "step" : step
            }
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
            "report_file" : report_file,
            "body_file" : body_file 
        }
    }

    # Create a unique input file name
    input_file = study_directory+ "/input_" + str(point_index)  + "_" + str(xyz_index) + "_" + f'{step: .2e}' + f'{cp_offset:.2e}' + ".json"
    write_input_file(input_dict, input_file)

    report = run_machline(input_file, run=RERUN_MACHLINE)

    # Pull norm of force
    CF = np.zeros(3)
    CF[0] = report["total_forces"]["Cx"]
    CF[1] = report["total_forces"]["Cy"]
    CF[2] = report["total_forces"]["Cz"]

    # if calc_adjoint, pull out norms (call them CF for convenience)
    if calc_adjoint:
        CF = np.zeros(3)
        CF[0] = report["norms_of_CF_sensitivities"]["norm_of_d_CFx"]
        CF[1] = report["norms_of_CF_sensitivities"]["norm_of_d_CFy"]
        CF[2] = report["norms_of_CF_sensitivities"]["norm_of_d_CFz"]


    # delete report
    os.remove(report_file)

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
    if delete_input:
        os.remove(input_filename)

    return report


# Function to generate dash styles based on the number of step sizes
def generate_dash_styles(num_steps):
    max_dash = 10  # Max dash length
    min_dash = 1   # Min dash length
    dash_step = (max_dash - min_dash) / (num_steps - 1)
    # Create dashes starting from max_dash to min_dash
    dash_styles = [(max_dash - i * dash_step, 2) for i in range(num_steps)]
    return dash_styles



def run_up_task(i, cp_offset, study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, mesh_name):
    CF_up = run_machline_for_cp_offset(cp_offset, study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, mesh_name)
    return CF_up

# Task to run "down" perturbations
def run_down_task(i, cp_offsets, study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, mesh_name):
    CF_down = run_machline_for_cp_offset(cp_offset, study_directory, calc_adjoint, perturb_point, point_index, xyz_index, -step, formulation, mesh_name)
    return CF_down


def k_loop(k, j, step, study_directory, calc_adjoint, cp_offsets, perturb_point, formulation, mesh_name):

    run_count_local = 0

    point_index = k
    xyz_index = j

    CF = run_machline_for_cp_offset(cp_offset, study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, mesh_name)
    run_count_local += 1

    # get up data
    CFx_up = CF[0]
    CFy_up = CF[1]
    CFz_up = CF[2]


    CF = run_machline_for_cp_offset(cp_offset, study_directory, calc_adjoint, perturb_point, point_index, xyz_index, -step, formulation, mesh_name)
    run_count_local += 1

    # get down data
    CFx_down = CF[0]
    CFy_down = CF[1]
    CFz_down = CF[2]

    # calculate sensitivity
    d_CFx = (CFx_up - CFx_down) / (2. * step) 
    d_CFy = (CFy_up - CFy_down) / (2. * step) 
    d_CFz = (CFz_up - CFz_down) / (2. * step) 
    
    return (k, j, d_CFx, d_CFy, d_CFz, run_count_local)


def process_in_batches(num_mesh_points, num_workers=60, j=1):
    run_count = 0
    results = []  # Store results for norms calculation
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for k in range(1, num_mesh_points + 1):
            # Submit jobs to the executor for the given axis (j)
            futures.append(executor.submit(k_loop, k, j, step, study_directory, calc_adjoint, cp_offset, perturb_point, formulation, mesh_name))

        # Collect results as they are completed
        for future in concurrent.futures.as_completed(futures):
            result = future.result()  # This will block until the future is complete
            results.append(result)  # Store the result
            k, j, d_CFx_result, d_CFy_result, d_CFz_result, run_count_local = result
            run_count += run_count_local

            
            d_CFx[(k + (j-1)*num_mesh_points) -1] = d_CFx_result  
            d_CFy[(k + (j-1)*num_mesh_points) -1] = d_CFy_result  
            d_CFz[(k + (j-1)*num_mesh_points) -1] = d_CFz_result  
            print(f"Completed mesh point {k} for axis {j}, run count: {run_count_local}")

    return d_CFx, d_CFy, d_CFz  # Return the accumulated results

    


if __name__=="__main__":

    # declare varaibles and inputs
    tStart = time.time()

    #####################################################################################
    #####################################################################################
    #### THINGS TO CHANGE FORA A NEW RUN #####

    mesh_name = "test_11"
    sonic = "super"
    num_mesh_points = 1190

    clones = [1,3,74,110,146,182,218,254,290,326,362,398,434,470,506,542,578,614,650,686,722,758,794,830,866,902,938,974,1010,1046,1082,1118,1154]
    # clones = [1]

    adjoint_cp_study = True

    cp_offset = 1.0e-5
    step = 1.0e-4 

    ###################################################################################
    ###################################################################################
    formulation = "dirichlet-source-free"


    study_directory = "studies/adjoint_studies/one_central_diff"
    
    # counter for number of times machline is run
    run_count = 0
    
    # final sensitivity vectors for each num_case cp_offset
    d_CFx = np.zeros(num_mesh_points*3)
    d_CFy = np.zeros(num_mesh_points*3)
    d_CFz = np.zeros(num_mesh_points*3)

    # print(cp_offsets)
    # sys.exit()

    original_directory = os.getcwd()

    makefile_directory = 'C:/Users/nathan/git-repos/MachLine'  # Adjust this path as needed

    # Change the working directory to the location of the Makefile
    os.chdir(makefile_directory)

    try:
        sp.run(['make', 'serial'])
        print("Machline make serial successful")
    except sp.CalledProcessError as e:
        print(f"Error during make: {e}")

    os.chdir(original_directory)

######################################### Central Diff Starts Here #############################################
    # perturb_point
    calc_adjoint = False
    perturb_point = True

    norms_data = []

    # xyz loop - iterate over axes outside of the process_in_batches function
    for j in range(1, 4):  # Loop over xyz axes


        d_CFx, d_CFy, d_CFz = process_in_batches(num_mesh_points, num_workers=60, j=j)
        print("size of d_CFx", len(d_CFx))

        
    
    
    d_CFx_norm = np.sqrt(np.sum(d_CFx * d_CFx))
    d_CFy_norm = np.sqrt(np.sum(d_CFy * d_CFy))
    d_CFz_norm = np.sqrt(np.sum(d_CFz * d_CFz))

    processed_d_CFx = []
    processed_d_CFy = []
    processed_d_CFz = []

    cloned_d_CFx = []
    cloned_d_CFy = []
    cloned_d_CFz = []


    for j in range(num_mesh_points):
        x_index = j
        print("x = ", x_index)
        y_index = num_mesh_points + j
        z_index = 2 * num_mesh_points + j
        processed_d_CFx.append((d_CFx[x_index], d_CFx[y_index], d_CFx[z_index]))
        processed_d_CFy.append((d_CFy[x_index], d_CFy[y_index], d_CFy[z_index]))
        processed_d_CFz.append((d_CFz[x_index], d_CFz[y_index], d_CFz[z_index]))

        excel_file = "studies/adjoint_studies/one_central_diff/excel_files/"+sonic+"_"+mesh_name + "_cp_"+f'{cp_offset:.2e}' + "_step_1e-" + str(step) + ".xlsx"

        if os.path.exists(excel_file):
            os.remove(excel_file)

        # Create a DataFrame for vector components (without norms)
        data_dict = {
            'd_CFx_x': [item[0] for item in processed_d_CFx],
            'd_CFx_y': [item[1] for item in processed_d_CFx],
            'd_CFx_z': [item[2] for item in processed_d_CFx],
            'd_CFy_x': [item[0] for item in processed_d_CFy],
            'd_CFy_y': [item[1] for item in processed_d_CFy],
            'd_CFy_z': [item[2] for item in processed_d_CFy],
            'd_CFz_x': [item[0] for item in processed_d_CFz],
            'd_CFz_y': [item[1] for item in processed_d_CFz],
            'd_CFz_z': [item[2] for item in processed_d_CFz],
        }

        # Convert to DataFrame and write to Excel
        df_vectors = pd.DataFrame(data_dict)
        df_vectors.to_excel(excel_file, index=False)


        # Now add the norms to the same file
        df_vectors['CFx_norm'] = d_CFx_norm
        df_vectors['CFy_norm'] = d_CFy_norm
        df_vectors['CFz_norm'] = d_CFz_norm

        # Save updated DataFrame with norms back to the same Excel file
        df_vectors.to_excel(excel_file, index=False)

        # Check if the current index (1-based) is in clone indices
        if (j + 1) in clones:
            # Store cloned points
            cloned_d_CFx.append((d_CFx[x_index], d_CFx[y_index], d_CFx[z_index]))
            cloned_d_CFy.append((d_CFy[x_index], d_CFy[y_index], d_CFy[z_index]))
            cloned_d_CFz.append((d_CFz[x_index], d_CFz[y_index], d_CFz[z_index]))

    processed_d_CFx.extend(cloned_d_CFx)
    processed_d_CFy.extend(cloned_d_CFy)
    processed_d_CFz.extend(cloned_d_CFz)
    print("                                            length of processed = ", len(processed_d_CFx))

    new_data = [processed_d_CFx, processed_d_CFy, processed_d_CFz]


    # set_path to a vtk containg just the points
    just_points_vtk = "studies/adjoint_studies/one_central_diff/vtk_files/just_points_"+mesh_name+".vtk"
    vtk_lines = read_vtk_file(just_points_vtk)

    # write central diff sensitivities to a vtk file file
    new_vtk = "studies/adjoint_studies/one_central_diff/vtk_files/super_"+mesh_name+"_cp_"+f'{cp_offset:.2e}' + "_step_1e-" + str(step) + ".vtk"
    if os.path.exists(new_vtk):
            os.remove(new_vtk)
    updated_vtk_content = add_vector_data_to_vtk(vtk_lines, new_data)
    write_vtk_file(new_vtk, updated_vtk_content)


    print("Norms: ",d_CFx_norm, d_CFy_norm, d_CFz_norm)


    tEnd = time.time()
    print("Elapsed time is {0} seconds".format(tEnd-tStart))

    print("Machline ran " + str(run_count) + " times")
    
    plt.show()
    sys.exit()
    
    