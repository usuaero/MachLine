import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, get_order_of_convergence, write_input_file, cases, line_styles
import studies.case_running_functions as crf


RERUN_MACHLINE = False
crf.RERUN_MH_ONLY = False
study_dir = "studies/supersonic_double_wedge_wing/"
plot_dir = study_dir + "plots/convergence/"


def run_quad_for_mach_aoa_and_mesh(M, alpha, grid, MCA):
    """Runs a case quad for the given Mach number, angle of attack, and mesh density."""

    # Parameters
    half_angle = 5.0

    # Storage locations
    case_name = "M_{0}_aoa_{1}_{2}_deg_{3}_MCA_{4}".format(M, alpha, int(half_angle), grid, MCA)
    mesh_file = study_dir + "meshes/diamond_{0}_deg_full_{1}.stl".format(int(half_angle), grid)
    results_file = study_dir + "results/"+case_name+".vtk"
    report_file = study_dir + "reports/"+case_name+".json"

    # Write out input file

    # Declare MachLine input
    input_dict = {
        "flow": {
            "freestream_velocity": [np.cos(np.radians(alpha)), 0.0, np.sin(np.radians(alpha))],
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "max_continuity_angle" : MCA,
            "force_sigma_match" : True,
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
    input_file = study_dir + "input.json"
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

    # Study parameters
    grids = ["coarse", "medium", "fine", "ultra_fine"]
    #grids = ["ultra_fine"]
    Ms = [1.5, 2.0, 3.0, 5.0]
    alphas = np.linspace(0.0, 5.0, 6)
    MCAs = [1.0, 45.0]

    # Loop through parameters
    N_sys = np.zeros(len(grids))
    l_avg = np.zeros(len(grids))
    C_F = np.zeros((len(grids), len(Ms), len(alphas), len(MCAs), 4, 3))
    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):
            for k, alpha in enumerate(alphas):
                for l, MCA in enumerate(MCAs):

                    N_sys[i], l_avg[i], C_F[i,j,k,l] = run_quad_for_mach_aoa_and_mesh(M, alpha, grid, MCA)

    # Set up plotting

    # Get errors
    err = np.abs((C_F[:-1] - C_F[-1])/C_F[-1])

    # We want the lowest MCA here
    l = 0
    max_continuity_angle = MCAs[l]

    # Plot C_x errors
    for j, M in enumerate(Ms):
        for k, alpha in enumerate(alphas):

            # Plot
            plt.figure()
            print()
            print("C_x")
            print("M={0}, alpha={1}".format(M, alpha))
            for m, case in enumerate(cases):
                plt.plot(l_avg[:-1], err[:,j,k,l,m,0], line_styles[m], label=case)

                # Get convergence
                print("Order for case {0}: {1}".format(case, get_order_of_convergence(l_avg, C_F[:,j,k,l,m,0], truth_from_results=True)))

            # Format
            plt.xlabel('$l_{avg}$')
            plt.ylabel('Fractional Error in $C_x$')
            plt.xscale('log')
            plt.yscale('log')
            #plt.title('$C_x$ error for $M={0},\,\\alpha={1}$'.format(M, alpha))
            if M==1.5 and alpha==1.0:
                plt.legend()
            plt.savefig(plot_dir+"err_C_x_M_{0}_alpha_{1}_MCA_{2}.pdf".format(M, alpha, max_continuity_angle))
            plt.savefig(plot_dir+"err_C_x_M_{0}_alpha_{1}_MCA_{2}.svg".format(M, alpha, max_continuity_angle))
            plt.close()

    # Plot C_z errors
    for j, M in enumerate(Ms):
        for k, alpha in enumerate(alphas):
            if k==0:
                continue

            # Plot
            plt.figure()
            print()
            print("C_z")
            print("M={0}, alpha={1}".format(M, alpha))
            for m, case in enumerate(cases):
                plt.plot(l_avg[:-1], err[:,j,k,l,m,2], line_styles[m], label=case)

                # Get convergence
                print("Order for case {0}: {1}".format(case, get_order_of_convergence(l_avg, C_F[:,j,k,l,m,2], truth_from_results=True)))

            # Format
            plt.xlabel('$l_{avg}$')
            plt.ylabel('Fractional Error in $C_z$')
            plt.xscale('log')
            plt.yscale('log')
            #plt.title('$C_z$ error for $M={0},\,\\alpha={1}$'.format(M, alpha))
            if M==1.5 and alpha==1.0:
                plt.legend()
            plt.savefig(plot_dir+"err_C_z_M_{0}_alpha_{1}_MCA_{2}.pdf".format(M, alpha, max_continuity_angle))
            plt.savefig(plot_dir+"err_C_z_M_{0}_alpha_{1}_MCA_{2}.svg".format(M, alpha, max_continuity_angle))
            plt.close()

    #Analyze convergence of Cx
    print()
    print("Cx Convergence Rate")
    print("-------------------")
    for i, MCA in enumerate(MCAs):
        print("---Max Continuity Angle: {0} degrees---".format(MCA))
        for l, case in enumerate(['ML', 'MH', 'SL', 'SH']):
            slopes = []
            for j, M in enumerate(Ms):
                for k, alpha in enumerate(alphas):
                    slopes.append(get_order_of_convergence(l_avg, C_F[:,j,k,i,l,0], truth_from_results=True))

            print("{0}: {1} +/- {2}".format(case, round(np.average(slopes), 4), round(np.std(slopes), 4)))

    #Analyze convergence of Cz
    print()
    print("Cz Convergence Rate")
    print("-------------------")
    for i, MCA in enumerate(MCAs):
        print("---Max Continuity Angle: {0} degrees---".format(MCA))
        for l, case in enumerate(['ML', 'MH', 'SL', 'SH']):
            slopes = []
            for j, M in enumerate(Ms):
                for k, alpha in enumerate(alphas):
                    if k==0:
                        continue
                    slopes.append(get_order_of_convergence(l_avg, C_F[:,j,k,i,l,2], truth_from_results=True))

            print("{0}: {1} +/- {2}".format(case, round(np.average(slopes), 4), round(np.std(slopes), 4)))

    # Morino-lower only plots

    # Plot C_x errors
    plt.figure()
    for j, M in enumerate(Ms):
        for k, alpha in enumerate(alphas):

            # Plot
            print()
            print("C_x")
            print("M={0}, alpha={1}".format(M, alpha))
            plt.plot(l_avg[:-1], err[:,j,k,0,0,0], 'k-')

            # Get convergence
            print("Order for case ML: {0}".format(get_order_of_convergence(l_avg, C_F[:,j,k,0,0,0], truth_from_results=True)))

    # Format
    plt.xlabel('$l_{avg}$')
    plt.ylabel('Fractional Error in $C_x$')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(plot_dir+"collected_err_C_x_ML.pdf".format(M, alpha))
    plt.savefig(plot_dir+"collected_err_C_x_ML.svg".format(M, alpha))
    plt.close()

    # Plot C_z errors
    plt.figure()
    for j, M in enumerate(Ms):
        for k, alpha in enumerate(alphas):
            if k==0:
                continue

            # Plot
            print()
            print("C_z")
            print("M={0}, alpha={1}".format(M, alpha))
            plt.plot(l_avg[:-1], err[:,j,k,0,0,2], 'k-')

            # Get convergence
            print("Order for case ML: {0}".format(get_order_of_convergence(l_avg, C_F[:,j,k,0,0,2], truth_from_results=True)))

    # Format
    plt.xlabel('$l_{avg}$')
    plt.ylabel('Fractional Error in $C_z$')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(plot_dir+"collected_err_C_z_ML.pdf".format(M, alpha))
    plt.savefig(plot_dir+"collected_err_C_z_ML.svg".format(M, alpha))
    plt.close()