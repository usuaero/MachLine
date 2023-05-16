import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, write_input_file, cases, line_styles
from studies.paraview_functions import extract_all_data, get_data_column_from_array, extract_plane_slice


RERUN_MACHLINE = False
study_dir = "studies/supersonic_double_wedge_wing/"
plot_dir = study_dir + "plots/shock_expansion/"
MCA = 1.0
half_angle = 5.0
gamma = 1.4


def run_machline(M, alpha, grid):
    # Runs MachLine with the given parameters

    # Storage locations
    case_name = "M_{0}_aoa_{1}_{2}_deg_{3}_MCA_{4}".format(M, alpha, int(half_angle), grid, MCA)
    mesh_file = study_dir + "meshes/diamond_{0}_deg_full_{1}.stl".format(int(half_angle), grid)
    results_file = study_dir + "results/"+case_name+".vtk"
    report_file = study_dir + "reports/"+case_name+".json"

    # Declare MachLine input
    input_dict = {
        "flow": {
            "freestream_velocity": [np.cos(np.radians(alpha)), 0.0, np.sin(np.radians(alpha))],
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

    return reports


def plot_comparison(M, alphas, grid, se_pressures, reports):
    # Plots the comparison of the diamond wing to shock-expansion theory

    # Loop through cases
    for i, case in enumerate(cases):

        case_name = "M_{0}_{1}_{2}".format(M, grid, case)

        # Parse out shock-expansion data
        Cp2 = se_pressures[:,0]
        Cp3 = se_pressures[:,1]
        Cp4 = se_pressures[:,2]
        Cp5 = se_pressures[:,3]

        # Initialize plot with SE data
        plt.figure()
        plt.plot(alphas, Cp2, "k-" )#, label="$C_{P_2}$")
        plt.plot(alphas, Cp3, "k--")#, label="$C_{P_3}$")
        plt.plot(alphas, Cp4, "k-.")#, label="$C_{P_4}$")
        plt.plot(alphas, Cp5, "k:" )#, label="$C_{P_5}$")

        # Initialize MachLine data storage
        C_p_ise_2_ML = np.zeros(len(alphas))
        C_p_ise_3_ML = np.zeros(len(alphas))
        C_p_ise_4_ML = np.zeros(len(alphas))
        C_p_ise_5_ML = np.zeros(len(alphas))
        C_p_2nd_2_ML = np.zeros(len(alphas))
        C_p_2nd_3_ML = np.zeros(len(alphas))
        C_p_2nd_4_ML = np.zeros(len(alphas))
        C_p_2nd_5_ML = np.zeros(len(alphas))
        C_p_sln_2_ML = np.zeros(len(alphas))
        C_p_sln_3_ML = np.zeros(len(alphas))
        C_p_sln_4_ML = np.zeros(len(alphas))
        C_p_sln_5_ML = np.zeros(len(alphas))
        C_p_lin_2_ML = np.zeros(len(alphas))
        C_p_lin_3_ML = np.zeros(len(alphas))
        C_p_lin_4_ML = np.zeros(len(alphas))
        C_p_lin_5_ML = np.zeros(len(alphas))
    
        # Loop through angles of attack to get data from MachLine
        for j, alpha in enumerate(alphas):

            # Get data
            result_file = reports[j][i]["input"]["output"]["body_file"]
            headers, data = extract_plane_slice(result_file, np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 0.0]), which_data='cell')
            x = get_data_column_from_array(headers, data, 'centroid:0')
            z = get_data_column_from_array(headers, data, 'centroid:2')
            C_p_ise = get_data_column_from_array(headers, data, 'C_p_ise')
            C_p_2nd = get_data_column_from_array(headers, data, 'C_p_2nd')
            C_p_sln = get_data_column_from_array(headers, data, 'C_p_sln')
            C_p_lin = get_data_column_from_array(headers, data, 'C_p_lin')

            # Figure out which points belong to which surfaces
            ind2 = np.where(np.logical_and(x < 0.5, z<0))
            ind3 = np.where(np.logical_and(x < 0.5, z>0))
            ind4 = np.where(np.logical_and(x > 0.5, z<0))
            ind5 = np.where(np.logical_and(x > 0.5, z>0))

            # Get average pressures
            C_p_ise_2_ML[j] = np.average(C_p_ise[ind2])
            C_p_ise_3_ML[j] = np.average(C_p_ise[ind3])
            C_p_ise_4_ML[j] = np.average(C_p_ise[ind4])
            C_p_ise_5_ML[j] = np.average(C_p_ise[ind5])
            C_p_2nd_2_ML[j] = np.average(C_p_2nd[ind2])
            C_p_2nd_3_ML[j] = np.average(C_p_2nd[ind3])
            C_p_2nd_4_ML[j] = np.average(C_p_2nd[ind4])
            C_p_2nd_5_ML[j] = np.average(C_p_2nd[ind5])
            C_p_sln_2_ML[j] = np.average(C_p_sln[ind2])
            C_p_sln_3_ML[j] = np.average(C_p_sln[ind3])
            C_p_sln_4_ML[j] = np.average(C_p_sln[ind4])
            C_p_sln_5_ML[j] = np.average(C_p_sln[ind5])
            C_p_lin_2_ML[j] = np.average(C_p_lin[ind2])
            C_p_lin_3_ML[j] = np.average(C_p_lin[ind3])
            C_p_lin_4_ML[j] = np.average(C_p_lin[ind4])
            C_p_lin_5_ML[j] = np.average(C_p_lin[ind5])

        # Plot data from MachLine
        plt.plot(alphas, C_p_ise_2_ML, 'kv', markersize=3, label="Isentropic")
        plt.plot(alphas, C_p_ise_3_ML, 'kv', markersize=3)
        plt.plot(alphas, C_p_ise_4_ML, 'kv', markersize=3)
        plt.plot(alphas, C_p_ise_5_ML, 'kv', markersize=3)
        plt.plot(alphas, C_p_2nd_2_ML, 'ks', markersize=3, label="Second-Order")
        plt.plot(alphas, C_p_2nd_3_ML, 'ks', markersize=3)
        plt.plot(alphas, C_p_2nd_4_ML, 'ks', markersize=3)
        plt.plot(alphas, C_p_2nd_5_ML, 'ks', markersize=3)
        plt.plot(alphas, C_p_sln_2_ML, 'ko', markersize=3, label="Slender-Body")
        plt.plot(alphas, C_p_sln_3_ML, 'ko', markersize=3)
        plt.plot(alphas, C_p_sln_4_ML, 'ko', markersize=3)
        plt.plot(alphas, C_p_sln_5_ML, 'ko', markersize=3)
        plt.plot(alphas, C_p_lin_2_ML, 'k^', markersize=3, label="Linear")
        plt.plot(alphas, C_p_lin_3_ML, 'k^', markersize=3)
        plt.plot(alphas, C_p_lin_4_ML, 'k^', markersize=3)
        plt.plot(alphas, C_p_lin_5_ML, 'k^', markersize=3)

        # Format
        plt.xlabel('$\\alpha$')
        plt.ylabel('$C_P$')
        plt.legend(fontsize=6, title_fontsize=6)
        plt.gca().invert_yaxis()
        plt.savefig(plot_dir+case_name+".pdf")
        plt.savefig(plot_dir+case_name+".svg")
        plt.close()


if __name__=="__main__":

    # options
    grids = ["coarse", "medium", "fine", "ultra_fine"]
    Ms = [1.5, 2.0, 3.0, 5.0]
    alphas = np.linspace(0.0, 5.0, 6)

    # Load shock-expansion data
    se_data = np.genfromtxt(study_dir + "shock_expansion_data.csv", delimiter=',')
    se_data = se_data.reshape((len(Ms), len(alphas), 4))

    # Loop through cases
    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):
            reports = []
            for k, alpha in enumerate(alphas):

                # Run MachLine
                reports.append(run_machline(M, alpha, grid))

            # Plot
            plot_comparison(M, alphas, grid, se_data[j,:,:], reports)