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
# from tqdm import tqdm 



RERUN_MACHLINE = True


def run_machline_for_cp_offset(cp_offset,study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, mesh_name, alpha=0, wake_present=False):
    
    print("Running MachLine for cp_offset = {0}".format(cp_offset))
    # Storage locations
    cp_offset_string = str(cp_offset).replace(".","_") + "_" + formulation.replace("-","_")
    case_name = "cp_offset_{0}".format(cp_offset_string)
    mesh_file = study_directory+"/meshes/"+mesh_name+".stl"
    results_file = study_directory+"/results/"+case_name+".vtk"
    # wake_file = study_directory+"/results/"+case_name+"_wake.vtk"
    report_file = study_directory+"/reports/1.json"
    control_point_file = study_directory+"/results/"+case_name+"_control_points.vtk"

    if (calc_adjoint):
        write_body_file = "body_file"
        body_file = study_directory + "/vtk_files/adjoint_"+mesh_name+"_cp_"+f'{cp_offset:.2e}'+".vtk"
    else:
        write_body_file = "xbody_file"
        body_file = "none"
    

    # create input file
    input_dict = {
        "flow": {
            "freestream_velocity": [1.0, 0.0, 0.1],
            "freestream_mach_number" : 0.5
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
            "report_file" : report_file,
            write_body_file : body_file 
        }
    }
    # Dump
    input_file = study_directory+"/input.json"
    write_input_file(input_dict, input_file)

    # Run MachLine
    report = run_machline(input_file, run=RERUN_MACHLINE)

    # # Pull out forces
    # C_F = np.zeros(3)
    # C_F[0] = report["total_forces"]["Cx"]
    # C_F[1] = report["total_forces"]["Cy"]
    # C_F[2] = report["total_forces"]["Cz"]

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


# task to be done in parallel
# def run_up_down_task(i,cp_offsets, study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, mesh_name):
    
#     # run "up" peruturbations
#     CF_up = run_machline_for_cp_offset(cp_offsets[i], study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, mesh_name)

#     # run "down" perturbations
#     CF_down = run_machline_for_cp_offset(cp_offsets[i], study_directory, calc_adjoint, perturb_point, point_index, xyz_index, -step, formulation, mesh_name)

#     return CF_up, CF_down



def k_loop(k, j, num_cp_offsets, step, study_directory, calc_adjoint, cp_offsets, perturb_point, formulation, mesh_name):
    CF = [None]*num_cp_offsets
    CFx_up = [None]*num_cp_offsets
    CFy_up = [None]*num_cp_offsets
    CFz_up = [None]*num_cp_offsets
    CFx_down = [None]*num_cp_offsets
    CFy_down = [None]*num_cp_offsets
    CFz_down = [None]*num_cp_offsets

    run_count_local = 0

    point_index = k
    xyz_index = j

    # perturb up
    for i in range(num_cp_offsets):
        CF[i] = run_machline_for_cp_offset(cp_offsets[i], study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, mesh_name)
        run_count_local += 1

    # get up data
    for i in range(num_cp_offsets):
        CFx_up[i] = CF[i][0]
        CFy_up[i] = CF[i][1]
        CFz_up[i] = CF[i][2]


     # perturb down
    for i in range(num_cp_offsets):
        CF[i] = run_machline_for_cp_offset(cp_offsets[i], study_directory, calc_adjoint, perturb_point, point_index, xyz_index, -step, formulation, mesh_name)
        run_count_local += 1

    # get down data
    for i in range(num_cp_offsets):
        CFx_down[i] = CF[i][0]
        CFy_down[i] = CF[i][1]
        CFz_down[i] = CF[i][2]

    # calculate sensitivity
    d_CFx = [(CFx_up[i] - CFx_down[i]) / (2. * step) for i in range(num_cp_offsets)]
    d_CFy = [(CFy_up[i] - CFy_down[i]) / (2. * step) for i in range(num_cp_offsets)]
    d_CFz = [(CFz_up[i] - CFz_down[i]) / (2. * step) for i in range(num_cp_offsets)]
    
    return (k, j, d_CFx, d_CFy, d_CFz, run_count_local)


def process_in_batches(num_mesh_points, num_workers=10):
    run_count = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for j in range(1, 4):  # Loop over xyz axes
            for k in range(1, num_mesh_points + 1):
                # Submit jobs to the executor
                futures.append(executor.submit(k_loop, k, j, num_cp_offsets, step, study_directory, calc_adjoint, cp_offsets, perturb_point, formulation, mesh_name))

        # Collect results as they are completed
        for future in concurrent.futures.as_completed(futures):
            result = future.result()  # This will block until the future is complete
            k, j, d_CFx, d_CFy, d_CFz, run_count_local = result
            run_count += run_count_local
            print(f"Completed mesh point {k} for axis {j}, run count: {run_count_local}")



if __name__=="__main__":

    # declare varaibles and inputs
    tStart = time.time()
    alpha = 0

    ##########################################
    #### THINGS TO CHANGE FORA A NEW RUN #####
    mesh_name = "test_11"
    sonic = "supersonic"
    num_mesh_points = 1190
    num_cp_offsets = 10
    step = 1.0e-4   # initial step size (gets smaller)
    initial_step_exp = 4
    num_step_size_runs = 7
    adjoint_cp_study = False
    
    #   - mesh_file 
    #   - just_points_vtk 
    #   - clones
    #   - num_mesh_points
    #   - step size, initial exp

    ########################################
    # clones 
    clones = [1,3,74,110,146,182,218,254,290,326,362,398,434,470,506,542,578,614,650,686,722,758,794,830,866,902,938,974,1010,1046,1082,1118,1154]
    # clones = [1]

    # get spread of cp offsets
    cp_offsets = np.logspace(-11,-1, num_cp_offsets+1)
    cp_offsets = cp_offsets[:-1]
    # cp_offsets = [1.0e-4] #, 7.5e-12, 1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6]

    # Function to generate dash styles based on the number of step sizes
    def generate_dash_styles(num_steps):
        max_dash = 10  # Max dash length
        min_dash = 1   # Min dash length
        dash_step = (max_dash - min_dash) / (num_steps - 1)
        # Create dashes starting from max_dash to min_dash
        dash_styles = [(max_dash - i * dash_step, 2) for i in range(num_steps)]
        return dash_styles

    # Generate dash styles based on the number of step sizes
    dash_styles = generate_dash_styles(num_step_size_runs)

    # Generate colors from gray scale (lightening as step size decreases)
    colors = [(0.3 + 0.4 * (i / (num_step_size_runs - 1)),) * 3 for i in range(num_step_size_runs)]


    # set_path to a vtk containg just the points
    just_points_vtk = "studies/adjoint_studies/FD_norm_cp_offset_parallel_subsonic/vtk_files/just_points_"+mesh_name+".vtk"
    vtk_lines = read_vtk_file(just_points_vtk)

    
    wake_present = True
    # wake_type = "panel"
    # wake_type = "filaments"
    # formulation = "neumann-mass-flux-VCP"
    formulation = "dirichlet-source-free"
    calc_adjoint = False

    study_directory = "studies/adjoint_studies/FD_norm_cp_offset_parallel_subsonic"
    CF = list(range(num_cp_offsets))
    CFx_up = list(range(num_cp_offsets))
    CFy_up = list(range(num_cp_offsets))
    CFz_up = list(range(num_cp_offsets))
    CFx_down = list(range(num_cp_offsets))
    CFy_down = list(range(num_cp_offsets))
    CFz_down = list(range(num_cp_offsets))
    d_CFx = list(range(num_cp_offsets))
    d_CFy = list(range(num_cp_offsets)) 
    d_CFz = list(range(num_cp_offsets))
    d_CF_norm = list(range(num_cp_offsets))
    d_CFx_norm = list(range(num_cp_offsets))
    d_CFy_norm = list(range(num_cp_offsets))
    d_CFz_norm = list(range(num_cp_offsets))
    d_CFx_norm_adjoint = list(range(num_cp_offsets))
    d_CFy_norm_adjoint = list(range(num_cp_offsets))
    d_CFz_norm_adjoint = list(range(num_cp_offsets))
    
    
    # counter for number of times machline is run
    run_count = 0
    
    # final sensitivity vectors for each num_case cp_offset
    for i in range(num_cp_offsets):
        d_CFx[i] = np.zeros(num_mesh_points*3)
        d_CFy[i] = np.zeros(num_mesh_points*3)
        d_CFz[i] = np.zeros(num_mesh_points*3)

    # print(cp_offsets)
    # sys.exit()

    # Set Times New Roman as the default font
    plt.rcParams["font.family"] = "Times New Roman"

    # init plots for values of d_CFx, d_CFy, d_CFz
    figx, ax1 = plt.subplots(figsize=(10,6))
    figy, ax2 = plt.subplots(figsize=(10,6))
    figz, ax3 = plt.subplots(figsize=(10,6))

    # do an adjoint run
    perturb_point = False
    point_index = 0
    xyz_index = 0

    if (adjoint_cp_study):
        calc_adjoint = True

        for i in range(num_cp_offsets):
            d_CF_norm[i] = run_machline_for_cp_offset(cp_offsets[i],study_directory, calc_adjoint, perturb_point, point_index, xyz_index, -step, formulation, mesh_name)
            run_count += 1

        # get data
        for i in range(num_cp_offsets):
            d_CFx_norm_adjoint[i] = d_CF_norm[i][0]
            d_CFy_norm_adjoint[i] = d_CF_norm[i][1]
            d_CFz_norm_adjoint[i] = d_CF_norm[i][2]
    
        
        
        # Plot for norms_d_CFx (adjoint)
        ax1.plot(cp_offsets, d_CFx_norm_adjoint, linestyle='-', color='black', label= "Adjoint")
        
        # Plot for norms_d_CFy (adjoint)
        ax2.plot(cp_offsets, d_CFy_norm_adjoint, linestyle='-', color='black', label= "Adjoint")
        
        # Plot for norms_d_CFz (adjoint)
        ax3.plot(cp_offsets, d_CFz_norm_adjoint, linestyle='-', color='black', label= "Adjoint")
    else:
        # set to zero
        for i in range(num_cp_offsets):
            d_CFx_norm_adjoint[i] = 0.0
            d_CFy_norm_adjoint[i] = 0.0
            d_CFz_norm_adjoint[i] = 0.0

    # perturb_point
    calc_adjoint = False
    perturb_point = True

    norms_data = []

    # step size loop
    for m in range (num_step_size_runs):

        if m == 0:
            step = step
        else:
            step = step/10

        # xyz loop
        for j in range(1,4):
            
            # Call process_in_batches 
            process_in_batches(num_mesh_points, 20)

            # # parallelize k loop
            # with concurrent.futures.ProcessPoolExecutor() as executor:
            #     futures = [executor.submit(k_loop, k, j, num_cp_offsets, step, study_directory, calc_adjoint, cp_offsets, perturb_point, formulation, mesh_name)
            #                for k in range(1, num_mesh_points + 1)]

            #     # collect results
            #     for future in concurrent.futures.as_completed(futures):
            #         k, j, d_CFx_result, d_CFy_result, d_CFz_result, run_count_local = future.result()

            #         run_count += run_count_local

            #         # store results
            #         for i in range(num_cp_offsets):
            #             d_CFx[i][(k + (j-1)*num_mesh_points) - 1] = d_CFx_result[i]
            #             d_CFy[i][(k + (j-1)*num_mesh_points) - 1] = d_CFy_result[i]
            #             d_CFz[i][(k + (j-1)*num_mesh_points) - 1] = d_CFz_result[i]
            
            # # 
            # vertex loop
            # for k in range(1, num_mesh_points + 1):

            #     point_index = k
            #     xyz_index = j

            #     # perturb up
            #     with ProcessPoolExecutor() as executor:
            #         # submit tasks for each cp offset
            #         futures = [executor.submit(run_up_down_task, i, cp_offsets, study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, mesh_name)
            #                    for i in range(num_cp_offsets)]
                    
            #         # collect data
            #         for i, future in enumerate(futures):

            #     # for i in range(num_cp_offsets):
            #     #     CF[i] = run_machline_for_cp_offset(cp_offsets[i],study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, mesh_name)
            #     #     run_count += 1


            #     # get up data
            #     for i in range(num_cp_offsets):
            #         CFx_up[i] = CF[i][0]
            #         CFy_up[i] = CF[i][1]
            #         CFz_up[i] = CF[i][2]


            #     # perturb down
            #     for i in range(num_cp_offsets):
            #         CF[i] = run_machline_for_cp_offset(cp_offsets[i],study_directory, calc_adjoint, perturb_point, point_index, xyz_index, -step, formulation, mesh_name)
            #         run_count += 1


            #     # get up data
            #     for i in range(num_cp_offsets):
            #         CFx_down[i] = CF[i][0]
            #         CFy_down[i] = CF[i][1]
            #         CFz_down[i] = CF[i][2]

            #     # calculate sensitivity wrt this design variable and store
            #     for i in range(num_cp_offsets):
            #         d_CFx[i][(k + (j-1)*num_mesh_points) -1] = (CFx_up[i] - CFx_down[i]) / (2.*step)  
            #         d_CFy[i][(k + (j-1)*num_mesh_points) -1] = (CFy_up[i] - CFy_down[i]) / (2.*step)  
            #         d_CFz[i][(k + (j-1)*num_mesh_points) -1] = (CFz_up[i] - CFz_down[i]) / (2.*step)  

            
        # for each step size, do the following: 

        # calc norms
        for i in range(num_cp_offsets):
            d_CFx_norm[i] = np.sqrt(np.sum(d_CFx[i][:] * d_CFx[i][:]))
            d_CFy_norm[i] = np.sqrt(np.sum(d_CFy[i][:] * d_CFy[i][:]))
            d_CFz_norm[i] = np.sqrt(np.sum(d_CFz[i][:] * d_CFz[i][:]))

            processed_d_CFx = []
            processed_d_CFy = []
            processed_d_CFz = []

            cloned_d_CFx = []
            cloned_d_CFy = []
            cloned_d_CFz = []

            print("len ", len(d_CFx))

            for j in range(num_mesh_points):
                x_index = j
                print("x = ", x_index)
                y_index = num_mesh_points + j
                z_index = 2 * num_mesh_points + j
                processed_d_CFx.append((d_CFx[i][x_index], d_CFx[i][y_index], d_CFx[i][z_index]))
                processed_d_CFy.append((d_CFy[i][x_index], d_CFy[i][y_index], d_CFy[i][z_index]))
                processed_d_CFz.append((d_CFz[i][x_index], d_CFz[i][y_index], d_CFz[i][z_index]))

                excel_file = "studies/adjoint_studies/FD_norm_cp_offset_parallel_subsonic/excel_files/subsonic_FD_"+mesh_name + "_cp_"+f'{cp_offsets[i]:.2e}' + "_step_1e-" + str(initial_step_exp+ m) + ".xlsx"

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
                df_vectors['CFx_norm'] = d_CFx_norm[i]
                df_vectors['CFy_norm'] = d_CFy_norm[i]
                df_vectors['CFz_norm'] = d_CFz_norm[i]

                # Save updated DataFrame with norms back to the same Excel file
                df_vectors.to_excel(excel_file, index=False)

                # Check if the current index (1-based) is in clone indices
                if (j + 1) in clones:
                    # Store cloned points
                    cloned_d_CFx.append((d_CFx[i][x_index], d_CFx[i][y_index], d_CFx[i][z_index]))
                    cloned_d_CFy.append((d_CFy[i][x_index], d_CFy[i][y_index], d_CFy[i][z_index]))
                    cloned_d_CFz.append((d_CFz[i][x_index], d_CFz[i][y_index], d_CFz[i][z_index]))

            processed_d_CFx.extend(cloned_d_CFx)
            processed_d_CFy.extend(cloned_d_CFy)
            processed_d_CFz.extend(cloned_d_CFz)

            new_data = [processed_d_CFx, processed_d_CFy, processed_d_CFz]

            # processed_d_CFx = np.array(processed_d_CFx)
            # processed_d_CFy = np.array(processed_d_CFy)
            # processed_d_CFz = np.array(processed_d_CFz)

            # write central diff sensitivities to a vtk file file
            new_vtk = "studies/adjoint_studies/FD_norm_cp_offset_parallel_subsonic/vtk_files/subsonic_FD_"+mesh_name+"_cp_"+f'{cp_offsets[i]:.2e}' + "_step_1e-" + str(initial_step_exp+ m) + ".vtk"
            if os.path.exists(new_vtk):
                    os.remove(new_vtk)
            updated_vtk_content = add_vector_data_to_vtk(vtk_lines, new_data)
            write_vtk_file(new_vtk, updated_vtk_content)


        
        # Plot for norms_d_CFx
        ax1.plot(cp_offsets, d_CFx_norm, linestyle=(0, dash_styles[m]), color="black", label= "Step size = 10^-" + str(initial_step_exp + m))
        
        # Plot for norms_d_CFy
        ax2.plot(cp_offsets, d_CFy_norm, linestyle=(0, dash_styles[m]), color="black", label= "Step size = 10^-" + str(initial_step_exp + m))
        
        # Plot for norms_d_CFz
        ax3.plot(cp_offsets, d_CFz_norm, linestyle=(0, dash_styles[m]), color="black", label= "Step size = 10^-" + str(initial_step_exp + m))
        

        # save data in result excel
        # Append data for this step to the norms_data list
        norms_data.append({
            "Step Size": f"10^-{initial_step_exp + m}",
            "CP Offset": cp_offsets,
            "d_CFx Norm": d_CFx_norm,
            "d_CFy Norm": d_CFy_norm,
            "d_CFz Norm": d_CFz_norm,
        })
    
    # finish CFx figure
    ax1.set_xscale("log")
    ax1.set_xlabel("Control Point Offset", fontsize = 14)
    ax1.set_ylabel("Norm of Sensitivity Vector", fontsize = 14)
    # ax1.set_title("Norm of CFx Sensitivities vs CP Offset ("+mesh_name + ", " +sonic+")", fontsize = 14)
    ax1.set_title("Norm of CFx Sensitivities vs CP Offset ("+mesh_name+", " +sonic+")", fontsize = 14)
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 12)
    # ax1.set_ylim(0.0, 5.0)
    ax1.set_yscale("log")
    ax1.tick_params(axis='both', labelsize=12)  # For CFx
    figx.subplots_adjust(right = 0.8)


    fig_file_fx = "/figures/subsonic_FD_"+mesh_name+"_"+str(num_cp_offsets) + "_cp_offsets_dCFx.pdf"
    figx.savefig(study_directory + fig_file_fx)

    # finish CFy figure
    ax2.set_xscale("log")
    ax2.set_xlabel("Control Point Offset", fontsize = 14)
    ax2.set_ylabel("Norm of Sensitivity Vector", fontsize = 14)
    # ax2.set_title("Norm of CFy Sensitivities vs CP Offset ("+mesh_name + ", " +sonic+")", fontsize = 14)
    ax2.set_title("Norm of CFy Sensitivities vs CP Offset ("+mesh_name+", " +sonic+")", fontsize = 14)
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 12)
    # ax2.set_ylim(0.0, 5.0)
    ax2.set_yscale("log")
    ax2.tick_params(axis='both', labelsize=12)  # For CFy
    figy.subplots_adjust
    fig_file_fy = "/figures/subsonic_FD_"+mesh_name+"_"+str(num_cp_offsets) + "_cp_offsets_dCFy.pdf"
    figy.savefig(study_directory + fig_file_fy)

    # finish CFy figure
    ax3.set_xscale("log")
    ax3.set_xlabel("Control Point Offset", fontsize = 14)
    ax3.set_ylabel("Norm of Sensitivity Vector", fontsize = 14)
    # ax3.set_title("Norm of CFz Sensitivities vs CP Offset ("+mesh_name + ", " +sonic+")", fontsize = 14)
    ax3.set_title("Norm of CFz Sensitivities vs CP Offset ("+mesh_name+", " +sonic+")", fontsize = 14)
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 12)
    # ax3.set_ylim(0.0, 5.0)
    ax3.set_yscale("log")
    ax3.tick_params(axis='both', labelsize=12)  # For CFz
    figz.subplots_adjust(right = 0.8)

    fig_file_fz = "/figures/subsonic_FD_"+mesh_name+"_"+str(num_cp_offsets) + "_cp_offsets_dCFz.pdf"
    figz.savefig(study_directory + fig_file_fz)

    with open("studies/adjoint_studies/FD_norm_cp_offset_parallel_subsonic/pickle_jar/subsonic_FD_"+mesh_name+"_"+str(num_cp_offsets) + "_cp_offsets_dCFx.pkl", 'wb') as f:
        pickle.dump(figx, f)


    with open("studies/adjoint_studies/FD_norm_cp_offset_parallel_subsonic/pickle_jar/subsonic_FD_"+mesh_name+"_"+str(num_cp_offsets) + "_cp_offsets_dCFy.pkl", 'wb') as f:
        pickle.dump(figy, f)


    with open("studies/adjoint_studies/FD_norm_cp_offset_parallel_subsonic/pickle_jar/subsonic_FD_"+mesh_name+"_"+str(num_cp_offsets) + "_cp_offsets_dCFz.pkl", 'wb') as f:
        pickle.dump(figz, f)


    # loading in the pickle
    #####################################################
    # Load figure from the file
    # with open("studies/adjoint_studies/fd_norm_cp_offset_parallel_subsonic/pickle_jar/central_diff_"+mesh_name+"_"+str(num_cp_offsets) + "_cp_offsets_dCFz.pkl", 'rb') as f:
    #     loaded_fig = pickle.load(f)

    # #Access the axes and modify the title
    # ax = loaded_fig.gca()  # Get current axes
    # ax.set_title("Silly Pickles")

    # # Display the modified figure
    # loaded_fig.show()
    #####################################################
        
    
    # write norms results to excel file in results
    # Create a DataFrame from the norms_data
    df_norms = pd.DataFrame(norms_data)

    # Convert the list of dictionaries to a suitable format
    # Here we expand the DataFrame for CP Offset and norms
    expanded_data = {
        "Step Size": [],
        "CP Offset": [],
        "d_CFx Norm": [],
        "d_CFy Norm": [],
        "d_CFz Norm": [],
        "d_CFx Norm Adjoint": [],
        "d_CFy Norm Adjoint": [],
        "d_CFz Norm Adjoint": []
    }

    # Populate the expanded_data for central difference norms
    for data in norms_data:
        for offset, cf_x, cf_y, cf_z in zip(data["CP Offset"], data["d_CFx Norm"], data["d_CFy Norm"], data["d_CFz Norm"]):
            expanded_data["Step Size"].append(data["Step Size"])
            expanded_data["CP Offset"].append(offset)
            expanded_data["d_CFx Norm"].append(cf_x)
            expanded_data["d_CFy Norm"].append(cf_y)
            expanded_data["d_CFz Norm"].append(cf_z)
            expanded_data["d_CFx Norm Adjoint"].append(None)  # Add adjoint norms
            expanded_data["d_CFy Norm Adjoint"].append(None)
            expanded_data["d_CFz Norm Adjoint"].append(None)

    # Now add the adjoint data without step size
    adjoint_size = len(d_CFx_norm_adjoint)
    offset_size = len(cp_offsets)
    print("adjoint_size = ", adjoint_size)
    print("cp_offset size = ", offset_size)
    print("adjoint", d_CFx_norm_adjoint[0])
    for offset, adj_cf_x, adj_cf_y, adj_cf_z in zip(cp_offsets, d_CFx_norm_adjoint, d_CFy_norm_adjoint, d_CFz_norm_adjoint):
        expanded_data["Step Size"].append("Adjoint")  # Mark the row as "Adjoint" since no step size applies
        expanded_data["CP Offset"].append(offset)
        expanded_data["d_CFx Norm"].append(None)  # No central difference data in this row
        expanded_data["d_CFy Norm"].append(None)  # Fill central difference norms with None
        expanded_data["d_CFz Norm"].append(None)
        expanded_data["d_CFx Norm Adjoint"].append(adj_cf_x)  # Add adjoint norms
        expanded_data["d_CFy Norm Adjoint"].append(adj_cf_y)
        expanded_data["d_CFz Norm Adjoint"].append(adj_cf_z)

    # Convert to DataFrame
    df_expanded = pd.DataFrame(expanded_data)

    # Save to Excel file
    excel_file = "studies/adjoint_studies/FD_norm_cp_offset_parallel_subsonic/results/subsonic_FD" + mesh_name + "_" + str(num_cp_offsets) + "_cp_offsets_norms_vs_cp_offset_data.xlsx"
    df_expanded.to_excel(excel_file, index=False)


    tEnd = time.time()
    print("Elapsed time is {0} seconds".format(tEnd-tStart))

    print("Machline ran " + str(run_count) + " times")
    
    plt.show()
    sys.exit()
    
    