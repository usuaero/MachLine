import os
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp


RERUN_MACHLINE = True


def run_machline_for_trefftz(trefftz):
    # Storage locations
    case_name = "NACA_trefftz_{0}".format(trefftz)
    mesh_file = "studies/supersonic_double_wedge_wing_filaments/meshes/naca_0010_AR_10_full_medium.stl"
    results_file = "studies/supersonic_double_wedge_wing_filaments/results/"+case_name+".vtk"
    wake_file = "studies/supersonic_double_wedge_wing_filaments/results/"+case_name+"_wake.vtk"
    report_file = "studies/supersonic_double_wedge_wing_filaments/reports/"+case_name+".json"

    # create input file
    input_dict = {
        "flow": {
            "freestream_velocity": [np.cos(np.radians(5)), 0.0, np.sin(np.radians(5))],
            "freestream_mach_number" : 0.25
            
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "reference" : {
                "area" : 10.0
            },
            "wake_model" : {
                "wake_present" : True,
                "append_wake" : True,
                "trefftz_distance": trefftz,
                "wake_type": "panel"
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
            "wake_file" : wake_file,
            "report_file" : report_file
        }
    }
    # Dump
    input_file = "studies/supersonic_double_wedge_wing_filaments/input.json"
    write_input_file(input_dict, input_file)

    # Run quad
    report = run_machline(input_file, run=RERUN_MACHLINE)

    # Pull out forces
    C_F = np.zeros(3)
    C_F[0] = report["total_forces"]["Cx"]
    C_F[1] = report["total_forces"]["Cy"]
    C_F[2] = report["total_forces"]["Cz"]

    # Get system dimension and average characteristic length
    N_sys = report["solver_results"]["system_dimension"]
    l_avg = report["mesh_info"]["average_characteristic_length"]

    return N_sys, l_avg, C_F

    



def get_json(file_path):
    # import json file from file path
    json_string = open(file_path).read()

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
    num_cases = 30
    N_sys = list(range(num_cases))
    l_avg =list(range(num_cases))
    C_F = list(range(num_cases))
    C_x = list(range(num_cases))
    C_y = list(range(num_cases))
    C_z = list(range(num_cases))

    # Run cases
    trefftz_distances = np.linspace(1,500,num_cases)
    for i in range(num_cases):
        N_sys[i], l_avg[i], C_F[i] = run_machline_for_trefftz(trefftz_distances[i])


    # get data
    for i in range(num_cases):
        C_x[i] = C_F[i][0]
        C_y[i] = C_F[i][1]
        C_z[i] = C_F[i][2]
        
    #  plot
    plt.figure()
    plt.plot(trefftz_distances,C_x, label="C_x")
    plt.plot(trefftz_distances,C_y,label="C_y")
    plt.plot(trefftz_distances,C_z,label="C_z")
    
    plt.xlabel("Trefftz Distance")
    plt.ylabel("$C_F$")
    plt.legend()
    plt.show()

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