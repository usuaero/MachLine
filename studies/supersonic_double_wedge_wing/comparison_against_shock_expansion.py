import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, write_input_file, cases, line_styles
from studies.paraview_functions import extract_all_data, get_data_column_from_array, extract_plane_slice


RERUN_MACHLINE = True
study_dir = "studies/supersonic_double_wedge_wing/"
plot_dir = study_dir + "plots/shock_expansion/"


def run_comparison(M, alpha, grid, MCA, se_pressures):
    # Runs the comparison of the diamond wing to shock-expansion theory

    # Parameters
    half_angle = 5.0
    R_G = 287.058
    gamma = 1.4
    T_inf = 300.0
    c_inf = np.sqrt(gamma*R_G*T_inf)
    rho = 1.225
    p_inf = 1.0e5

    # Storage locations
    case_name = "M_{0}_aoa_{1}_{2}_deg_{3}_MCA_{4}".format(M, alpha, int(half_angle), grid, MCA)
    mesh_file = study_dir + "meshes/diamond_{0}_deg_full_{1}.stl".format(int(half_angle), grid)
    results_file = study_dir + "results/"+case_name+".vtk"
    report_file = study_dir + "reports/"+case_name+".json"

    # Declare MachLine input
    input_dict = {
        "flow": {
            "freestream_velocity": [M*c_inf*np.cos(np.radians(alpha)), 0.0, M*c_inf*np.sin(np.radians(alpha))],
            "gamma" : gamma,
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "max_continuity_angle" : MCA,
            "wake_model": {
                "append_wake" : False,
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
                "second-order" : True,
                "isentropic" : True,
                "slender-body" : True,
                "linear" : True
            }
        },
        "output" : {
            "body_file" : results_file,
            "report_file" : report_file
        }
    }

    # Dump
    input_file = study_dir + "diamond_input.json"
    write_input_file(input_dict, input_file)

    # Run
    reports = run_quad(input_file, run=RERUN_MACHLINE)
    
    # Parse out shock-expansion data
    Cp2 = se_pressures[0]
    Cp3 = se_pressures[1]
    Cp4 = se_pressures[2]
    Cp5 = se_pressures[3]
    
    # Loop through reports
    for i, (case, report) in enumerate(zip(cases, reports)):

        # Get data
        result_file = report["input"]["output"]["body_file"]
        headers, data = extract_plane_slice(result_file, np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 0.0]), which_data='cell')
        x = get_data_column_from_array(headers, data, 'centroid:0')
        C_p_ise = get_data_column_from_array(headers, data, 'C_p_ise')
        C_p_2nd = get_data_column_from_array(headers, data, 'C_p_2nd')
        C_p_sln = get_data_column_from_array(headers, data, 'C_p_sln')
        C_p_lin = get_data_column_from_array(headers, data, 'C_p_lin')

        # Plot data from MachLine
        plt.figure()
        plt.plot(x, C_p_2nd, 'ks', markersize=3, label='Second-Order')
        plt.plot(x, C_p_ise, 'kv', markersize=3, label='Isentropic')
        plt.plot(x, C_p_sln, 'ko', markersize=3, label='Slender-Body')
        plt.plot(x, C_p_lin, 'k^', markersize=3, label='Linear')

        # Plot data from shock-expansion theory
        x = np.linspace(0.0, 1.0, 100)
        Cp_upper = np.ones_like(x)
        Cp_upper[:50] *= Cp3
        Cp_upper[50:] *= Cp5
        Cp_lower = np.ones_like(x)
        Cp_lower[:50] *= Cp2
        Cp_lower[50:] *= Cp4
        plt.plot(x, Cp_upper, 'k-', label='Shock-Expansion', linewidth=0.5)
        plt.plot(x, Cp_lower, 'k-', linewidth=0.5)

        # Format
        plt.xlabel('$x$')
        plt.ylabel('$C_P$')
        plt.legend(fontsize=6, title_fontsize=6)
        plt.savefig(plot_dir+case_name+"_{0}.pdf".format(case))
        plt.savefig(plot_dir+case_name+"_{0}.svg".format(case))
        plt.close()


if __name__=="__main__":

    # options
    grids = ["coarse", "medium", "fine", "ultra_fine"]
    Ms = [1.5, 2.0, 3.0, 5.0]
    alphas = np.linspace(0.0, 5.0, 6)
    MCAs = [1.0, 2.5, 45.0]

    # Load shock-expansion data
    se_data = np.genfromtxt(study_dir + "shock_expansion_data.csv", delimiter=',')
    se_data = se_data.reshape((len(grids), len(Ms), len(alphas), len(MCAs), 4))

    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):
            for k, alpha in enumerate(alphas):
                for l, MCA in enumerate(MCAs):

                    run_comparison(M, alpha, grid, MCA, se_data[i,j,k,l])