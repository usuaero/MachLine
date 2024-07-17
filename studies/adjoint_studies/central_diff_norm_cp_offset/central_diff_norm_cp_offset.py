import os
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp
import time
import sys



RERUN_MACHLINE = True


def run_machline_for_cp_offset(cp_offset,study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation, alpha=0, wake_present=False):
    
    print("Running MachLine for cp_offset = {0}".format(cp_offset))
    # Storage locations
    cp_offset_string = str(cp_offset).replace(".","_") + "_" + formulation.replace("-","_")
    case_name = "cp_offset_{0}".format(cp_offset_string)
    mesh_file = study_directory+"/meshes/sphere_coarse2.stl"
    results_file = study_directory+"/results/"+case_name+".vtk"
    # wake_file = study_directory+"/results/"+case_name+"_wake.vtk"
    report_file = study_directory+"/reports/"+case_name+".json"
    control_point_file = study_directory+"/results/"+case_name+"_control_points.vtk"
    

    # create input file
    input_dict = {
        "flow": {
            "freestream_velocity": [np.cos(np.radians(alpha)),np.sin(np.radians(alpha)), 0.0 ]
        },
        "geometry": {
            "file": mesh_file,
            "wake_model" : {
                "wake_present" : wake_present
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
                "incompressible" : True
            }
        },
        "output" : {
            # "body_file" : results_file,
            "report_file" : report_file,
            # "control_point_file" : control_point_file
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


    # Get system dimension and average characteristic length
    N_sys = report["solver_results"]["system_dimension"]
    l_avg = report["mesh_info"]["average_characteristic_length"]

    # delete report
    # os.remove(report)  doesnt work

    return N_sys, l_avg, CF

    



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


if __name__=="__main__":

    # declare varaibles and inputs
    tStart = time.time()
    alpha = 0
    num_cases = 1000
    # mach = 2.0
    wake_present = False
    # wake_type = "panel"
    # wake_type = "filaments"
    # formulation = "neumann-mass-flux-VCP"
    formulation = "dirichlet-source-free"
    calc_adjoint = False

    study_directory = "studies/adjoint_studies/central_diff_norm_cp_offset"
    N_sys = list(range(num_cases))
    l_avg =list(range(num_cases))
    CF = list(range(num_cases))
    CFx_up = list(range(num_cases))
    CFy_up = list(range(num_cases))
    CFz_up = list(range(num_cases))
    CFx_down = list(range(num_cases))
    CFy_down = list(range(num_cases))
    CFz_down = list(range(num_cases))
    d_CFx = list(range(num_cases))
    d_CFy = list(range(num_cases)) 
    d_CFz = list(range(num_cases))
    d_CF_norm = list(range(num_cases))
    d_CFx_norm = list(range(num_cases))
    d_CFy_norm = list(range(num_cases))
    d_CFz_norm = list(range(num_cases))
    
    # counter for number of times machline is run
    run_count = 0

    # number of mesh points/vertices
    points = 26
    step = 1.0e-9
    
    # final sensitivity vectors for each num_case cp_offset
    for i in range(num_cases):
        d_CFx[i] = np.zeros(points*3)
        d_CFy[i] = np.zeros(points*3)
        d_CFz[i] = np.zeros(points*3)

    # get spread of cp offsets
    cp_offsets = np.logspace(-12,0, num_cases+1)
    cp_offsets = cp_offsets[:-1]
    print(cp_offsets)
    # sys.exit()

    # init plots for values of d_CFx, d_CFy, d_CFz
    figx, ax1 = plt.subplots(figsize=(10,6))
    figy, ax2 = plt.subplots(figsize=(10,6))
    figz, ax3 = plt.subplots(figsize=(10,6))

    # do an adjoint run
    calc_adjoint = True
    perturb_point = False
    point_index = 0
    xyz_index = 0

    for i in range(num_cases):
        N_sys[i], l_avg[i], d_CF_norm[i] = run_machline_for_cp_offset(cp_offsets[i],study_directory, calc_adjoint, perturb_point, point_index, xyz_index, -step, formulation)
        run_count += 1

    # get data
    for i in range(num_cases):
        d_CFx_norm[i] = d_CF_norm[i][0]
        d_CFy_norm[i] = d_CF_norm[i][1]
        d_CFz_norm[i] = d_CF_norm[i][2]

    # Plot for norms_d_CFx (adjoint)
    ax1.plot(cp_offsets, d_CFx_norm, linestyle='--', label= "Adjoint")
    
    # Plot for norms_d_CFy (adjoint)
    ax2.plot(cp_offsets, d_CFy_norm, linestyle='--', label= "Adjoint")
    
    # Plot for norms_d_CFz (adjoint)
    ax3.plot(cp_offsets, d_CFz_norm, linestyle='--', label= "Adjoint")
    

    # perturb_point
    calc_adjoint = False
    perturb_point = True

    # step size loop
    for m in range (9):

        if m == 0:
            step = step
        else:
            step = step*10

        # xyz loop
        for j in range(1,4):
            
            # vertex loop
            for k in range(1, points + 1):

                point_index = k
                xyz_index = j

                # perturb up
                for i in range(num_cases):
                    N_sys[i], l_avg[i], CF[i] = run_machline_for_cp_offset(cp_offsets[i],study_directory, calc_adjoint, perturb_point, point_index, xyz_index, step, formulation)
                    run_count += 1


                # get up data
                for i in range(num_cases):
                    CFx_up[i] = CF[i][0]
                    CFy_up[i] = CF[i][1]
                    CFz_up[i] = CF[i][2]


                # perturb down
                for i in range(num_cases):
                    N_sys[i], l_avg[i], CF[i] = run_machline_for_cp_offset(cp_offsets[i],study_directory, calc_adjoint, perturb_point, point_index, xyz_index, -step, formulation)
                    run_count += 1


                # get up data
                for i in range(num_cases):
                    CFx_down[i] = CF[i][0]
                    CFy_down[i] = CF[i][1]
                    CFz_down[i] = CF[i][2]

                # calculate sensitivity wrt this design variable and store
                for i in range(num_cases):
                    d_CFx[i][(k + (j-1)*points) -1] = (CFx_up[i] - CFx_down[i]) / (2.*step)  
                    d_CFy[i][(k + (j-1)*points) -1] = (CFy_up[i] - CFy_down[i]) / (2.*step)  
                    d_CFz[i][(k + (j-1)*points) -1] = (CFz_up[i] - CFz_down[i]) / (2.*step)  

            
        # for each central diff step size, do the following: 

        # calc norms
        for i in range(num_cases):
            d_CFx_norm[i] = np.sqrt(np.sum(d_CFx[i][:] * d_CFx[i][:]))
            d_CFy_norm[i] = np.sqrt(np.sum(d_CFy[i][:] * d_CFy[i][:]))
            d_CFz_norm[i] = np.sqrt(np.sum(d_CFz[i][:] * d_CFz[i][:]))
        
        # Plot for norms_d_CFx
        ax1.plot(cp_offsets, d_CFx_norm, label= "Step size = 1.0e-" + str(9 - m))
        
        # Plot for norms_d_CFy
        ax2.plot(cp_offsets, d_CFy_norm, label= "Step size = 1.0e-" + str(9 - m))
        
        # Plot for norms_d_CFz
        ax3.plot(cp_offsets, d_CFz_norm, label= "Step size = 1.0e-" + str(9 - m))
        
    
    
    # finish CFx figure
    ax1.set_xscale("log")
    ax1.set_xlabel("Control Point Offset")
    ax1.set_ylabel("Norm of Sensitivity Vector")
    ax1.set_title("Norm of CFx Sensitivities vs CP Offset (Central Diff)")
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_ylim(0.0, 5.0)
    figx.subplots_adjust(right = 0.8)

    fig_file_fx = "/figures/central_diff_norm_cp_offset_coarse_sphere_CFx_sensitivities.png"
    figx.savefig(study_directory + fig_file_fx)

    # finish CFy figure
    ax2.set_xscale("log")
    ax2.set_xlabel("Control Point Offset")
    ax2.set_ylabel("Norm of Sensitivity Vector")
    ax2.set_title("Norm of CFy Sensitivities vs CP Offset (Central Diff)")
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_ylim(0.0, 5.0)
    figy.subplots_adjust(right = 0.8)

    fig_file_fy = "/figures/central_diff_norm_cp_offset_coarse_sphere_CFy_sensitivities.png"
    figy.savefig(study_directory + fig_file_fy)

    # finish CFy figure
    ax3.set_xscale("log")
    ax3.set_xlabel("Control Point Offset")
    ax3.set_ylabel("Norm of Sensitivity Vector")
    ax3.set_title("Norm of CFz Sensitivities vs CP Offset (Central Diff)")
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax3.set_ylim(0.0, 5.0)
    figz.subplots_adjust(right = 0.8)

    fig_file_fz = "/figures/central_diff_norm_cp_offset_coarse_sphere_CFz_sensitivities.png"
    figz.savefig(study_directory + fig_file_fz)
    
    


    tEnd = time.time()
    print("Elapsed time is {0} seconds".format(tEnd-tStart))

    print("Machline ran " + str(run_count) + " times")
    
    plt.show()
    
    
    