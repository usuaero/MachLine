import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, write_input_file, cases, line_styles
from studies.paraview_functions import extract_all_data, get_data_column_from_array


RERUN_MACHLINE = True
study_dir = "studies/supersonic_cone/"
plot_dir = study_dir + "plots/"


def run_comparison(M, grid, half_angle):
    # Runs the comparison of a cone with the Taylor-MacColl result

    # Parameters
    gamma = 1.4

    # Storage locations
    case_name = "M_{0}_{1}_deg_{2}".format(M, int(half_angle), grid)
    #mesh_file = study_dir + "meshes/cone_{0}_deg_{1}_improved_full.vtk".format(int(half_angle), grid)
    mesh_file = study_dir + "meshes/cone_{0}_deg_{1}_improved.vtk".format(int(half_angle), grid)
    #mesh_file = study_dir + "meshes/cone_{0}_deg_{1}.vtk".format(int(half_angle), grid)
    body_file = study_dir + "results/"+case_name+".vtk"
    report_file = study_dir + "reports/"+case_name+".json"

    # Declare MachLine input
    input_dict = {
        "flow": {
            "freestream_velocity": [-1.0, 0.0, 0.0],
            "gamma" : gamma,
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+z",
            "mirror_about" : "xy",
            "max_continuity_angle" : 45.0,
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
            "body_file" : body_file,
            "report_file" : report_file
        }
    }

    # Dump
    input_file = study_dir + "cone_input.json"
    write_input_file(input_dict, input_file)

    # Run
    reports = run_quad(input_file, run=RERUN_MACHLINE)

    # Extract data
    C_p_2nd_avg = np.zeros(4)
    C_p_2nd_std = np.zeros(4)
    C_p_ise_avg = np.zeros(4)
    C_p_ise_std = np.zeros(4)
    C_p_lin_avg = np.zeros(4)
    C_p_lin_std = np.zeros(4)
    C_p_sln_avg = np.zeros(4)
    C_p_sln_std = np.zeros(4)
    t = np.zeros(4)
    for i, report in enumerate(reports):
    
        try:

            # Get data
            t[i] = report["total_runtime"]
            headers, data = extract_all_data(report["input"]["output"]["body_file"], which_data='cell')
            C_p_2nd = get_data_column_from_array(headers, data, 'C_p_2nd')
            C_p_ise = get_data_column_from_array(headers, data, 'C_p_ise')
            C_p_lin = get_data_column_from_array(headers, data, 'C_p_lin')
            C_p_sln = get_data_column_from_array(headers, data, 'C_p_sln')
            x = get_data_column_from_array(headers, data, 'centroid:0')

            # Calculate statistics
            C_p_2nd_avg[i] = np.average(C_p_2nd).item()
            C_p_2nd_std[i] = np.std(C_p_2nd).item()
            C_p_ise_avg[i] = np.average(C_p_ise).item()
            C_p_ise_std[i] = np.std(C_p_ise).item()
            C_p_sln_avg[i] = np.average(C_p_sln).item()
            C_p_sln_std[i] = np.std(C_p_sln).item()
            C_p_lin_avg[i] = np.average(C_p_lin).item()
            C_p_lin_std[i] = np.std(C_p_lin).item()
            
        except:
            t[i] = np.nan
            C_p_2nd_avg[i] = np.nan
            C_p_2nd_std[i] = np.nan
            C_p_ise_avg[i] = np.nan
            C_p_ise_std[i] = np.nan
            C_p_sln_avg[i] = np.nan
            C_p_sln_std[i] = np.nan
            C_p_lin_avg[i] = np.nan
            C_p_lin_std[i] = np.nan

    return C_p_2nd_avg, C_p_2nd_std, C_p_ise_avg, C_p_ise_std, C_p_sln_avg, C_p_sln_std, C_p_lin_avg, C_p_lin_std, t


def get_analytic_data(filename):
    """Reads in the analytic data from the given file."""

    # Get raw data
    raw_data = np.genfromtxt(filename, dtype=str, delimiter=',')

    # Parse data
    Ms = raw_data[2:,1].astype(float)
    thetas =  raw_data[1,2:].astype(float)
    Cps = raw_data[2:,2:].astype(float)

    return Ms, thetas, Cps


if __name__=="__main__":

    # Study parameters
    grids = ["coarse", "medium", "fine"]
    #grids = ["fine"]
    Ms = [1.4, 1.5, 1.7, 2.0, 2.4, 2.8, 3.3, 4.0]
    half_angles = [2.5, 5, 10, 15]

    # Initialize storage
    C_p_2nd_avg = np.zeros((len(grids), len(Ms), len(half_angles), 4))
    C_p_ise_avg = np.zeros((len(grids), len(Ms), len(half_angles), 4))
    C_p_sln_avg = np.zeros((len(grids), len(Ms), len(half_angles), 4))
    C_p_lin_avg = np.zeros((len(grids), len(Ms), len(half_angles), 4))
    C_p_2nd_std = np.zeros((len(grids), len(Ms), len(half_angles), 4))
    C_p_ise_std = np.zeros((len(grids), len(Ms), len(half_angles), 4))
    C_p_sln_std = np.zeros((len(grids), len(Ms), len(half_angles), 4))
    C_p_lin_std = np.zeros((len(grids), len(Ms), len(half_angles), 4))

    # Get analytic data
    Ms_anl, thetas_anl, Cps_anl = get_analytic_data(study_dir + "Cone Data Zero AoA.csv")

    # Run cases
    t = np.zeros((len(grids), len(Ms), len(half_angles), 4))
    for i, grid in enumerate(grids):
        for k, half_angle in enumerate(half_angles):

            # Get all Mach numbers
            for j, M in enumerate(Ms):

                if M == 4.0 and half_angle == 15.0:
                    C_p_2nd_avg[i,j,k,:] = np.nan
                    C_p_ise_avg[i,j,k,:] = np.nan
                    C_p_sln_avg[i,j,k,:] = np.nan
                    C_p_lin_avg[i,j,k,:] = np.nan
                    C_p_2nd_std[i,j,k,:] = np.nan
                    C_p_ise_std[i,j,k,:] = np.nan
                    C_p_sln_std[i,j,k,:] = np.nan
                    C_p_lin_std[i,j,k,:] = np.nan
                    t[i,j,k,:] = np.nan
                else:
                    C_p_2nd_avg[i,j,k], C_p_2nd_std[i,j,k], C_p_ise_avg[i,j,k], C_p_ise_std[i,j,k], C_p_sln_avg[i,j,k], C_p_sln_std[i,j,k], C_p_lin_avg[i,j,k], C_p_lin_std[i,j,k], t[i,j,k,:] = run_comparison(M, grid, half_angle)

            # Loop through case plots
            for j, (case, line_style) in enumerate(zip(cases, line_styles)):

                # Plot MachLine data
                plt.figure()
                plt.errorbar(Ms, C_p_2nd_avg[i,:,k,j], fmt='ks', yerr=C_p_2nd_std[i,:,k,j], elinewidth=0.5, markersize=3, label='Second-Order')
                plt.errorbar(Ms, C_p_ise_avg[i,:,k,j], fmt='kv', yerr=C_p_ise_std[i,:,k,j], elinewidth=0.5, markersize=3, label='Isentropic')
                plt.errorbar(Ms, C_p_sln_avg[i,:,k,j], fmt='ko', yerr=C_p_sln_std[i,:,k,j], elinewidth=0.5, markersize=3, label='Slender-Body')
                plt.errorbar(Ms, C_p_lin_avg[i,:,k,j], fmt='k^', yerr=C_p_lin_std[i,:,k,j], elinewidth=0.5, markersize=3, label='Linear')

                # Get analytic data for this half angle
                i_theta = np.where(thetas_anl == half_angle)
                plt.plot(Ms_anl[:7], Cps_anl[:7,i_theta].flatten(), 'k', label='Taylor-MacColl', linewidth=1)

                plt.xlabel("$M_\infty$")
                plt.ylabel("$C_P$")
                plt.ylim(bottom=0.0, top=1.1*np.nanmax(C_p_lin_avg[i,:,k,j]).item())
                if case=='MH':
                    plt.legend(fontsize=6, title_fontsize=6)
                plt.savefig(plot_dir + "C_p_over_M_{0}_deg_{1}_{2}.pdf".format(half_angle, case, grid))
                plt.savefig(plot_dir + "C_p_over_M_{0}_deg_{1}_{2}.svg".format(half_angle, case, grid))
                plt.close()

    # Print timing
    t_avg = np.nanmean(t, axis=1)
    for i, grid in enumerate(grids):
        for j, half_angle in enumerate(half_angles):
            for k, case in enumerate(cases):
                print(grid, half_angle, case, t_avg[i,j,k])

    ## Determine max standard deviation
    #std_max_ise = np.nanmax(np.nanmax(np.nanmax(C_p_ise_s_dev))).item()
    #std_max_2nd = np.nanmax(np.nanmax(np.nanmax(C_p_2nd_s_dev))).item()
    #std_max_lin = np.nanmax(np.nanmax(np.nanmax(C_p_lin_s_dev))).item()
    #std_max_sln = np.nanmax(np.nanmax(np.nanmax(C_p_sln_s_dev))).item()
    #print()
    #print("Maximum standard deviations:")
    #print("    Isentropic: ", std_max_ise)
    #print("    Second-order: ", std_max_2nd)
    #print("    Linear: ", std_max_lin)
    #print("    Slender-body: ", std_max_sln)

    # Plot convergence
    for k, half_angle in enumerate(half_angles):
        for j, M in enumerate(Ms):

            # Calculate error
            err = np.abs((C_p_ise_avg[:-1,j,k,:] - C_p_ise_avg[-1,j,k,:])/C_p_ise_avg[-1,j,k,:])

            plt.figure()
            for i, case in enumerate(cases):
                print(err[:,i])
                plt.plot(t[:-1,j,k,i], err[:,i], line_styles[i], label=case)
            plt.xlabel("Run Time $[s]$")
            plt.ylabel("Fractional Error in $C_{P_{ise}}$")
            plt.xscale('log')
            plt.yscale('log')
            plt.savefig(plot_dir + "accuracy_vs_time_{0}_M_{1}.pdf".format(int(half_angle), M))
            plt.savefig(plot_dir + "accuracy_vs_time_{0}_M_{1}.svg".format(int(half_angle), M))
            plt.close()

    ## Plot convergence
    #plt.figure()
    #for k, half_angle in enumerate(half_angles):
    #    if k==1:
    #        plt.plot(Ms, C_p_ise_avg[0,:,k], mfc='none', marker='o', ls='', mec='k', markersize=10, label='Coarse')
    #        plt.plot(Ms, C_p_ise_avg[1,:,k], mfc='none', marker='o', ls='', mec='k', markersize=7, label='Medium')
    #        plt.plot(Ms, C_p_ise_avg[2,:,k], mfc='none', marker='o', ls='', mec='k', markersize=4, label='Fine')
    #    else:
    #        plt.plot(Ms, C_p_ise_avg[0,:,k], mfc='none', marker='o', ls='', mec='k', markersize=10)
    #        plt.plot(Ms, C_p_ise_avg[1,:,k], mfc='none', marker='o', ls='', mec='k', markersize=7)
    #        plt.plot(Ms, C_p_ise_avg[2,:,k], mfc='none', marker='o', ls='', mec='k', markersize=4)
    #plt.xlabel("$M_\infty$")
    #plt.ylabel("$C_P$")
    #plt.ylim(bottom=0.0)
    #plt.legend(fontsize=6, title_fontsize=6)
    #plt.savefig("studies/supersonic_cone_flow_study/plots/C_p_ise_convergence.pdf")

    #plt.figure()
    #for k, half_angle in enumerate(half_angles):
    #    if k==1:
    #        plt.plot(Ms, C_p_2nd_avg[0,:,k], mfc='none', marker='o', ls='', mec='k', markersize=10, label='Coarse')
    #        plt.plot(Ms, C_p_2nd_avg[1,:,k], mfc='none', marker='o', ls='', mec='k', markersize=5, label='Medium')
    #        plt.plot(Ms, C_p_2nd_avg[2,:,k], mfc='none', marker='o', ls='', mec='k', markersize=3, label='Fine')
    #    else:
    #        plt.plot(Ms, C_p_2nd_avg[0,:,k], mfc='none', marker='o', ls='', mec='k', markersize=10)
    #        plt.plot(Ms, C_p_2nd_avg[1,:,k], mfc='none', marker='o', ls='', mec='k', markersize=5)
    #        plt.plot(Ms, C_p_2nd_avg[2,:,k], mfc='none', marker='o', ls='', mec='k', markersize=3)
    #plt.xlabel("$M_\infty$")
    #plt.ylabel("$C_P$")
    #plt.ylim(bottom=0.0)
    #plt.legend(fontsize=6, title_fontsize=6)
    #plt.savefig("studies/supersonic_cone_flow_study/plots/C_p_2nd_convergence.pdf")

    #plt.figure()
    #for k, half_angle in enumerate(half_angles):
    #    if k==1:
    #        plt.plot(Ms, C_p_sln_avg[0,:,k], mfc='none', marker='o', ls='', mec='k', markersize=10, label='Coarse')
    #        plt.plot(Ms, C_p_sln_avg[1,:,k], mfc='none', marker='o', ls='', mec='k', markersize=5, label='Medium')
    #        plt.plot(Ms, C_p_sln_avg[2,:,k], mfc='none', marker='o', ls='', mec='k', markersize=3, label='Fine')
    #    else:
    #        plt.plot(Ms, C_p_sln_avg[0,:,k], mfc='none', marker='o', ls='', mec='k', markersize=10)
    #        plt.plot(Ms, C_p_sln_avg[1,:,k], mfc='none', marker='o', ls='', mec='k', markersize=5)
    #        plt.plot(Ms, C_p_sln_avg[2,:,k], mfc='none', marker='o', ls='', mec='k', markersize=3)
    #plt.xlabel("$M_\infty$")
    #plt.ylabel("$C_P$")
    #plt.ylim(bottom=0.0)
    #plt.legend(fontsize=6, title_fontsize=6)
    #plt.savefig("studies/supersonic_cone_flow_study/plots/C_p_sln_convergence.pdf")

    #plt.figure()
    #for k, half_angle in enumerate(half_angles):
    #    if k==1:
    #        plt.plot(Ms, C_p_lin_avg[0,:,k], mfc='none', marker='o', ls='', mec='k', markersize=10, label='Coarse')
    #        plt.plot(Ms, C_p_lin_avg[1,:,k], mfc='none', marker='o', ls='', mec='k', markersize=5, label='Medium')
    #        plt.plot(Ms, C_p_lin_avg[2,:,k], mfc='none', marker='o', ls='', mec='k', markersize=3, label='Fine')
    #    else:
    #        plt.plot(Ms, C_p_lin_avg[0,:,k], mfc='none', marker='o', ls='', mec='k', markersize=10)
    #        plt.plot(Ms, C_p_lin_avg[1,:,k], mfc='none', marker='o', ls='', mec='k', markersize=5)
    #        plt.plot(Ms, C_p_lin_avg[2,:,k], mfc='none', marker='o', ls='', mec='k', markersize=3)
    #plt.xlabel("$M_\infty$")
    #plt.ylabel("$C_P$")
    #plt.ylim(bottom=0.0)
    #plt.legend(fontsize=6, title_fontsize=6)
    #plt.savefig("studies/supersonic_cone_flow_study/plots/C_p_lin_convergence.pdf")