import os
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp
import time



RERUN_MACHLINE = True


def run_machline_for_cp_offset(cp_offset,study_directory, alpha=0, wake_present=False, formulation="dirichlet-source-free"):
    
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
            }
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

    # Pull norm of force sensi
    d_CF = list(range(3))
    d_CF[0] = report["CF_sensitivities"]["d_CFx"]
    d_CF[1] = report["CF_sensitivities"]["d_CFy"]
    d_CF[2] = report["CF_sensitivities"]["d_CFz"]


    # Get system dimension and average characteristic length
    N_sys = report["solver_results"]["system_dimension"]
    l_avg = report["mesh_info"]["average_characteristic_length"]

    # delete report
    # os.remove(report)

    return N_sys, l_avg, d_CF

    



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
    calc_adjoint = True

    study_directory = "studies/adjoint_studies/adjoint_vals_cp_offset"
    N_sys = list(range(num_cases))
    l_avg =list(range(num_cases))
    d_CF = list(range(num_cases))
    d_CFx = list(range(num_cases))
    d_CFy = list(range(num_cases))
    d_CFz = list(range(num_cases))

    points = 26


    # Run cases
    cp_offsets = np.logspace(-12,0, num_cases)
    
    
    for i in range(num_cases):
        N_sys[i], l_avg[i], d_CF[i] = run_machline_for_cp_offset(cp_offsets[i],study_directory, alpha=alpha, wake_present=wake_present, formulation=formulation)


    # get data
    for i in range(num_cases):
        d_CFx[i] = d_CF[i][0]
        d_CFy[i] = d_CF[i][1]
        d_CFz[i] = d_CF[i][2]
    

    d_CFx_coord = np.zeros(num_cases)
    d_CFy_coord = np.zeros(num_cases)
    d_CFz_coord = np.zeros(num_cases)


    # Plot for values of d_CFx_dx5
    figx, ax1 = plt.subplots()
    figy, ax2 = plt.subplots()
    figz, ax3 = plt.subplots()

    for j in range(points):

        for k in range(3):

            for i in range(num_cases):
                d_CFx_coord[i] = d_CFx[i][j+(k+1)]
                d_CFy_coord[i] = d_CFy[i][j+(k+1)]
                d_CFz_coord[i] = d_CFz[i][j+(k+1)]
            

            if k == 0: 
                xyz_string = "x"
            elif k == 1: 
                xyz_string = "y"
            elif k == 2: 
                xyz_string = "z"
            # Plot for norms_d_CFx
            # ax1.plot(cp_offsets, d_CFx_coord, label= "wrt_" + xyz_string + str(j+1))
            ax1.plot(cp_offsets, d_CFx_coord)
            

            # Plot for norms_d_CFy
            # ax2.plot(cp_offsets, d_CFy_coord, label= "wrt_" + xyz_string + str(j+1))
            ax2.plot(cp_offsets, d_CFy_coord)
            

            # Plot for norms_d_CFz
            # ax3.plot(cp_offsets, d_CFz_coord, label= "wrt_" + xyz_string + str(j+1))
            ax3.plot(cp_offsets, d_CFz_coord)
            
    
    
    # finish CFx figure
    ax1.set_xscale("log")
    ax1.set_xlabel("Control Point Offset")
    ax1.set_ylabel("CFx Sensitivities")
    # ax1.legend()
    ax1.set_ylim(-1.0, 1.0)

    fig_file_fx = "/figures/adjoint_vals_cp_offset_coarse_sphere_CFx_sensitivities.png"
    figx.savefig(study_directory + fig_file_fx)
    figx.show()

    # finish CFy figure
    ax2.set_xscale("log")
    ax2.set_xlabel("Control Point Offset")
    ax2.set_ylabel("CFy Sensitivities")
    # ax2.legend()
    ax2.set_ylim(-1.0, 1.0)

    fig_file_fy = "/figures/adjoint_vals_cp_offset_coarse_sphere_CFy_sensitivities.png"
    figy.savefig(study_directory + fig_file_fy)
    figy.show()

    # finish CFy figure
    ax3.set_xscale("log")
    ax3.set_xlabel("Control Point Offset")
    ax3.set_ylabel("CFz Sensitivities")
    # ax3.legend()
    ax3.set_ylim(-1.0, 1.0)

    fig_file_fz = "/figures/adjoint_vals_cp_offset_coarse_sphere_CFz_sensitivities.png"
    figz.savefig(study_directory + fig_file_fz)
    figz.show()


    tEnd = time.time()
    print("Elapsed time is {0} seconds".format(tEnd-tStart))
    
    
    
    