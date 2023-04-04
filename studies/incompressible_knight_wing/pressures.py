import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, write_input_file


RERUN_MACHLINE = True


def run_quad_for_aoa_and_mesh(alpha, density):
    """Runs a case quad for the given angle of attack and mesh density."""

    # Storage locations
    case_name = "aoa_{0}_{1}".format(alpha, density)
    plot_dir = "studies/incompressible_knight_wing/plots/"
    mesh_file = "studies/incompressible_knight_wing/meshes/knight_wing_{0}.stl".format(density)
    results_file = "studies/incompressible_knight_wing/results/"+case_name+".vtk"
    report_file = "studies/incompressible_knight_wing/reports/"+case_name+".json"
    data_file = 'studies/incompressible_knight_wing/data/'+case_name+'.csv'

    # Write out input file

    # Declare MachLine input
    input_dict = {
        "flow": {
            "freestream_velocity": [np.cos(np.radians(alpha)), 0.0, np.sin(np.radians(alpha))]
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "mirror_about" : 'xz',
            "wake_model": {
            },
            "reference": {
                "area": 4.0
            }
        },
        "solver": {
            "formulation": "morino"
        },
        "post_processing" : {
            "pressure_rules" : {
                "incompressible" : True
            }
        },
        "output" : {
            "body_file" : results_file,
            "report_file" : report_file
        }
    }

    # Dump
    input_file = "studies/incompressible_knight_wing/input.json"
    write_input_file(input_dict, input_file)

    # Run quad
    reports = run_quad(input_file, run=RERUN_MACHLINE)

    # Pull out forces
    C_F = np.zeros((4,3))
    for i, report in enumerate(reports):
        C_F[i,0] = report["total_forces"]["Cx"]
        C_F[i,1] = report["total_forces"]["Cy"]
        C_F[i,2] = report["total_forces"]["Cz"]

    # Get system dimension and average characteristic length
    N_sys = reports[0]["solver_results"]["system_dimension"]
    l_avg = reports[0]["mesh_info"]["average_characteristic_length"]

    return N_sys, l_avg, C_F


if __name__=="__main__":

    # Run cases
    alphas = np.linspace(-8.0, 15.0 , 5)
    densities = ['coarse', 'medium', 'fine']
    for alpha in alphas:
        for density in densities:
            run_quad_for_aoa_and_mesh(alpha, density)

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