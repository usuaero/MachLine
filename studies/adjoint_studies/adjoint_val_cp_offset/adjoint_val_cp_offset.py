import os
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp
import time



RERUN_MACHLINE = True


def run_machline_for_cp_offset(cp_offset,study_directory, coord, alpha=0, wake_present=False, formulation="dirichlet-source-free"):
    
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

    # Pull norm of force sensitivities
    d_CF_d_coord = np.zeros(3)
    d_CF_d_coord[0] = report["CF_sensitivities"]["d_CFx"][coord]
    d_CF_d_coord[1] = report["CF_sensitivities"]["d_CFy"][coord]
    d_CF_d_coord[2] = report["CF_sensitivities"]["d_CFz"][coord]

    # Get system dimension and average characteristic length
    N_sys = report["solver_results"]["system_dimension"]
    l_avg = report["mesh_info"]["average_characteristic_length"]

    return N_sys, l_avg, d_CF_d_coord

    



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
    num_cases = 10
    # mach = 2.0
    wake_present = False
    # wake_type = "panel"
    # wake_type = "filaments"
    # formulation = "neumann-mass-flux-VCP"
    formulation = "dirichlet-source-free"
    calc_adjoint = True

    # sensitivity wrt which point?
    point = 10

    # senstivities wrt a point's x y or z coordinate?
    xyz = 1
    
    # the minus 1 converts to python indexing
    coord = xyz*point - 1

    if xyz == 1: 
        xyz_string = "x"
    elif xyz == 2: 
        xyz_string = "y"
    elif xyz == 3: 
        xyz_string = "z"



    study_directory = "studies/adjoint_studies/adjoint_val_cp_offset"
    N_sys = list(range(num_cases))
    l_avg =list(range(num_cases))
    d_CF_d_coord = list(range(num_cases))
    d_CFx_d_coord = list(range(num_cases))
    d_CFy_d_coord = list(range(num_cases))
    d_CFz_d_coord = list(range(num_cases))

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

    cp_offsets = np.logspace(-12,0, num_cases)
    
    
    for i in range(num_cases):
        N_sys[i], l_avg[i], d_CF_d_coord[i] = run_machline_for_cp_offset(cp_offsets[i],study_directory, coord, alpha=alpha, wake_present=wake_present, formulation=formulation)


    # get data
    for i in range(num_cases):
        d_CFx_d_coord[i] = d_CF_d_coord[i][0]
        d_CFy_d_coord[i] = d_CF_d_coord[i][1]
        d_CFz_d_coord[i] = d_CF_d_coord[i][2]
        
    tEnd = time.time()
    print("Elapsed time is {0} seconds".format(tEnd-tStart))
    
    # Plot for values of d_CFx_dx5
    plt.figure()
    plt.plot(cp_offsets, d_CFx_d_coord, color='black', label="d_CFx_wrt_" + xyz_string + str(point))
    plt.xscale("log")
    plt.xlabel("Control Point Offset")
    plt.ylabel("Sensitivity of CFx wrt " + xyz_string + " " + str(point))
    plt.legend()
    plt.ylim(-1.0, 1.0)

    fig_file_fx = "/figures/adjoint_val_cp_offset_coarse_sphere_d_CF_wrt_" + xyz_string + str(point) + ".png"
    plt.savefig(study_directory + fig_file_fx)
    plt.show()

    # Plot for norms_d_CFy
    plt.figure()
    plt.plot(cp_offsets, d_CFy_d_coord, color='black', label="d_CFy_wrt_" + xyz_string + str(point))
    plt.xscale("log")
    plt.xlabel("Control Point Offset")
    plt.ylabel("Sensitivity of CFy wrt " + xyz_string + " " + str(point))
    plt.legend()
    plt.ylim(-1.0, 1.0)

    fig_file_fy = "/figures/adjoint_val_cp_offset_coarse_sphere_d_CF_wrt_" + xyz_string + str(point) + ".png"
    plt.savefig(study_directory + fig_file_fy)
    plt.show()

    # Plot for norms_d_CFz
    plt.figure()
    plt.plot(cp_offsets, d_CFz_d_coord, color='black', label="d_CFz_wrt_" + xyz_string + str(point))
    plt.xscale("log")
    plt.xlabel("Control Point Offset")
    plt.ylabel("Sensitivity of CFz wrt " + xyz_string + " " + str(point))
    plt.legend()
    plt.ylim(-1.0, 1.0)

    fig_file_fz = "/figures/adjoint_val_cp_offset_coarse_sphere_d_CF_wrt_" + xyz_string + str(point) + ".png"
    plt.savefig(study_directory + fig_file_fz)
    plt.show()
    
    
    
    
    
    
    
    # #  plot
    # plt.figure()
    # plt.plot(cp_offsets,norms_d_CFx, label="norm of CFx sensitivities")
    # plt.plot(cp_offsets,norms_d_CFy,label="norm of CFy sensitivities")
    # plt.plot(cp_offsets,norms_d_CFz,label="norm of CFz sensitivities")
    # plt.xscale("log")
    # plt.xlabel("Control Point Offset")
    # plt.ylabel("Norm of sensitivities")
    # plt.legend()

    # # set axis limits
    # plt.ylim(0.0,5.0)    
    

    # # get figure file name
    # fig_file = "/adjoint_val_cp_offset - coarse sphere - {0}.png".format(formulation)
    # plt.savefig(study_directory+fig_file)
    # plt.show()



    ## Locations of pressure slices
    ## Station 1 : 14.75 in
    ## Station 2 : 14.25 in
    ## Station 3 : 13.5 in
    ## Station 4 : 12.5 in
    ## Station 5 : 11.0 in
    ## Station 6 : 9.0 in
    ## Station 7 : 6.0 in
    ## Station 8 : 3.0 in

    ## Read in experimental data
    #raw_data = np.genfromtxt("studies/incompressible_knight_wing/experimental pressures -8.csv", delimiter=',', skip_header=2)

    ## Get pressure corrections
    #C_p_refs = raw_data[0,33::2]

    ## Clean up data
    #C_p_lower = []
    #x_lower = []
    #C_p_upper = []
    #x_upper = []
    #for i in range(8):

    #    # Get pressure correction
    #    #C_p_ref = raw_data[0,2*i+31]

    #    # Get data
    #    x_upper.append(raw_data[:,4*i])
    #    C_p_upper.append(raw_data[:,4*i+1])
    #    x_lower.append(raw_data[:,4*i+2])
    #    C_p_lower.append(raw_data[:,4*i+3])

    #    # Remove nans and correct pressure
    #    x_upper[i] = x_upper[i][np.logical_not(np.isnan(x_upper[i]))]
    #    x_lower[i] = x_lower[i][np.logical_not(np.isnan(x_lower[i]))]
    #    C_p_upper[i] = C_p_upper[i][np.logical_not(np.isnan(C_p_upper[i]))] - C_p_refs[i]
    #    C_p_lower[i] = C_p_lower[i][np.logical_not(np.isnan(C_p_lower[i]))] - C_p_refs[i]

    #    # Load panel data
    #    panel_data = np.genfromtxt("station {0} morino.csv".format(i+1), delimiter=',', skip_header=1)
    #    C_p_panel = panel_data[:,0]
    #    x_panel = panel_data[:,1]
    #    x_min = min(x_panel)
    #    x_max = max(x_panel)
    #    x_panel /= x_max-x_min # Normalize
    #    x_panel *= -1.0 # Flip
    #    x_panel += 0.25 # Shift

    #    # plot
    #    plt.figure(figsize=(5.5, 5))
    #    plt.plot(x_upper[i], C_p_upper[i], 'sk', label="Exp. Upper")
    #    plt.plot(x_lower[i], C_p_lower[i], 'vk', label="Exp. Lower")
    #    plt.plot(x_panel, C_p_panel, 'k--', label="Panel Method")
    #    plt.xlabel("$x/c$")
    #    plt.ylabel("$C_p$")
    #    plt.legend()
    #    plt.gca().invert_yaxis()
    #    plt.savefig("knight -8 station {0}.pdf".format(i+1))
    #    plt.savefig("knight -8 station {0}.svg".format(i+1))
    #    plt.close()