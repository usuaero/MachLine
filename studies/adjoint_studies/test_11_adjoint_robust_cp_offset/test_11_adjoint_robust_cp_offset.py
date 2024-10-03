import os
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp
import time



RERUN_MACHLINE = True


def run_machline_for_cp_offset(cp_offset,study_directory,mach = 0.0,alpha=0.0, wake_present=True, formulation="dirichlet-source-free"):
    
    print("Running MachLine for cp_offset = {0}".format(cp_offset))
    # Storage locations
    cp_offset_string = str(cp_offset).replace(".","_") + "_" + formulation.replace("-","_")
    case_name = "cp_offset_{0}".format(cp_offset_string)
    mesh_file = study_directory+"/meshes/test_mesh_11.stl"
    results_file = study_directory+"/results/"+case_name+".vtk"
    wake_file = study_directory+"/results/"+case_name+"_wake.vtk"
    report_file = study_directory+"/reports/"+case_name+".json"
    control_point_file = study_directory+"/results/"+case_name+"_control_points.vtk"
    

    # create input file
    input_dict = {
        "flow": {
            "freestream_velocity": [np.cos(np.radians(alpha)),0.0, np.sin(np.radians(alpha))],
            "freestream_mach_number" : mach
        },
        "geometry": {
            "file": mesh_file,
            "wake_model" : {
                "wake_present" : True,
                "append_wake" : True,
                "xtrefftz_distance": 19.0
                 },
            "adjoint_sensitivities" : {
                "calc_adjoint" : True,
                "robust_cp_dir" : True
            }
        },
    
        "solver": {
            "formulation" : formulation,
            "control_point_offset": cp_offset,
            "matrix_solver" : "GMRES"
        },
        "post_processing" : {
            "pressure_rules" : {
                "incompressible" : False,
                "isentropic" : True
            }
        },
        "output" : {
            "verbose" : False,
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

    # Pull norm of force sensitivities
    norms_d_CF = np.zeros(3)
    norms_d_CF[0] = report["norms_of_CF_sensitivities"]["norm_of_d_CFx"]
    norms_d_CF[1] = report["norms_of_CF_sensitivities"]["norm_of_d_CFy"]
    norms_d_CF[2] = report["norms_of_CF_sensitivities"]["norm_of_d_CFz"]

    # Get system dimension and average characteristic length
    N_sys = report["solver_results"]["system_dimension"]
    l_avg = report["mesh_info"]["average_characteristic_length"]

    os.remove(report_file)

    return N_sys, l_avg, norms_d_CF

    



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
    alpha = 5.0
    num_cases = 20
    mach = 2.0
    wake_present = True
    # wake_type = "panel"
    # wake_type = "filaments"
    # formulation = "neumann-mass-flux-VCP"
    formulation = "dirichlet-source-free"
    calc_adjoint = True



    study_directory = "studies/adjoint_studies/test_11_adjoint_robust_cp_offset"
    N_sys = list(range(num_cases))
    l_avg =list(range(num_cases))
    norms_d_CF = list(range(num_cases))
    norms_d_CFx = list(range(num_cases))
    norms_d_CFy = list(range(num_cases))
    norms_d_CFz = list(range(num_cases))

    # Run cases

    ################################################################
    # powers_of_10 = np.logspace(-12,0,13)
    # times_5 = 5.*np.logspace(-12,0,13)
    # cp_offsets = np.sort(np.concatenate((powers_of_10, times_5)))
    
    # # get rid of any duplicate numbers
    # cp_offsets = np.unique(cp_offsets)

    # # remove the last number because I don't want it
    # cp_offsets = cp_offsets[:-1]
    ################################################################

    cp_offsets = np.logspace(-12,-1, num_cases)
    
    
    for i in range(num_cases):
        N_sys[i], l_avg[i], norms_d_CF[i] = run_machline_for_cp_offset(cp_offsets[i],study_directory, mach=mach, alpha=alpha, wake_present=wake_present, formulation=formulation)


    # get data
    for i in range(num_cases):
        norms_d_CFx[i] = norms_d_CF[i][0]
        norms_d_CFy[i] = norms_d_CF[i][1]
        norms_d_CFz[i] = norms_d_CF[i][2]
        
    tEnd = time.time()
    print("Elapsed time is {0} seconds".format(tEnd-tStart))
    
    # Plot for norms_d_CFx
    plt.figure()
    plt.plot(cp_offsets, norms_d_CFx, color='black', label="norm of CFx sensitivities")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Control Point Offset")
    plt.ylabel("Norm of CFx sensitivities")
    plt.legend()
    # plt.ylim(0.0, 5.0)

    fig_file_fx = "/figures/adjoint_norm_cp_offset_test_11_CFx.png"
    plt.savefig(study_directory + fig_file_fx)
    plt.show()

    # Plot for norms_d_CFy
    plt.figure()
    plt.plot(cp_offsets, norms_d_CFy, color='black', label="norm of CFy sensitivities")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Control Point Offset")
    plt.ylabel("Norm of CFy sensitivities")
    plt.legend()
    # plt.ylim(0.0, 5.0)

    fig_file_fy = "/figures/adjoint_norm_cp_offset_test_11_CFy.png"
    plt.savefig(study_directory + fig_file_fy)
    plt.show()

    # Plot for norms_d_CFz
    plt.figure()
    plt.plot(cp_offsets, norms_d_CFz, color='black', label="norm of CFz sensitivities")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Control Point Offset")
    plt.ylabel("Norm of CFz sensitivities")
    plt.legend()
    # plt.ylim(0.0, 5.0)

    fig_file_fz = "/figures/adjoint_norm_cp_offset_test_11_CFz.png"
    plt.savefig(study_directory + fig_file_fz)
    plt.show()
    
    
    
    
    
    