import numpy as np
import matplotlib.pyplot as plt


study_dir = "studies/matrix_solvers/"
iteration_dir = study_dir + "iterations/"
plot_dir = study_dir + "plots/"
data_dir = study_dir + "data/"


def get_iterative_method_residual_history(mesh_root_name, solver, refinement, preconditioner, sort_system):
    # Reads in the residual history from the file

    # Get file location
    iter_file = iteration_dir + "{0}{1}_{2}_{3}_{4}_prec_history.csv".format(mesh_root_name, solver, refinement, preconditioner, "sorted" if sort_system else "unsorted")

    # Load data
    iteration_data = np.genfromtxt(iter_file, skip_header=4, delimiter=',', dtype=float)
    if solver == "GMRES":
        residual_history = iteration_data[:,1]
    else:
        residual_history = iteration_data[:,2]

    return residual_history


def get_overall_info(mesh_root_name, solver, refinement, preconditioner, sort_system):
    # Gets the run time and residual info for the given case from file

    # Get file location
    data_file = data_dir + "{0}solver_test_data.csv".format(mesh_root_name)

    # Read in file
    with open(data_file, 'r') as data_handle:
        lines = data_handle.readlines()

    # Loop through lines to find the info we need
    runtimes = []
    residuals = []
    for line in lines[1:]:
        split_line = line.split(',')
        if split_line[0] == solver and split_line[1] == refinement and split_line[2] == preconditioner and split_line[3].title() == str(sort_system):
            runtimes.append(float(split_line[6]))
            residuals.append(float(split_line[7]))

    return np.average(runtimes).item(), np.average(residuals).item()


if __name__=="__main__":

    # Options
    case_names = ["cone_10_deg_", "diamond_5_deg_full_", "full_config_"]
    iterative_solver_options = ["BJAC", "BSSOR", "GMRES"]
    direct_solver_options = ["QRUP", "FQRUP"]
    l_styles = ['-', '--', '-.']
    direct_styles = ['ks', 'kv']
    refinement_options = ["coarse", "medium", "fine"]
    colors = ['black', 'dimgray', 'silver']
    preconditioner_options = ["DIAG", "none"]
    sort_system_options = [False, True]

    # Run through
    for case_name in case_names:
        for sort_system in sort_system_options:
            for k, preconditioner in enumerate(preconditioner_options):

                # Standard iteration plots

                # Start the figure for this case
                plt.figure()

                # Loop through refinement levels
                for j, refinement in enumerate(refinement_options):

                    # Get iterative solver iterations
                    for i, solver in enumerate(iterative_solver_options):

                        res = get_iterative_method_residual_history(case_name, solver, refinement, preconditioner, sort_system)
                        if j==0:
                            plt.plot(res, linestyle=l_styles[i], color=colors[j], label=solver)
                        else:
                            plt.plot(res, linestyle=l_styles[i], color=colors[j])

                plt.xlabel('Iteration')
                plt.ylabel('Norm of Residual')
                plt.yscale('log')
                plt.gca().set_ylim(top=10.0, bottom=1.0e-13)
                plt.legend()
                plt.savefig(plot_dir+"{0}{1}_{2}_iteration_history.pdf".format(case_name, "sorted" if sort_system else "unsorted", preconditioner))
                plt.savefig(plot_dir+"{0}{1}_{2}_iteration_history.svg".format(case_name, "sorted" if sort_system else "unsorted", preconditioner))
                plt.close()

                # Time vs accuracy plots

                # Loop through refinement levels
                for j, refinement in enumerate(refinement_options):

                    # Start figure
                    plt.figure()

                    # Get iterative solver iterations
                    for i, solver in enumerate(iterative_solver_options):

                        # Get run info
                        res = get_iterative_method_residual_history(case_name, solver, refinement, preconditioner, sort_system)
                        runtime, residual = get_overall_info(case_name, solver, refinement, preconditioner, sort_system)

                        # Parse iterations to run time
                        t = np.arange(1, len(res)+1, dtype=float)
                        t /= t[-1]
                        t *= runtime

                        # Plot
                        plt.plot(t, res, linestyle=l_styles[i], color='k', label=solver)

                    # Get direct solver results
                    for i, solver in enumerate(direct_solver_options):
                        runtime, residual = get_overall_info(case_name, solver, refinement, preconditioner, sort_system)
                        plt.plot(runtime, residual, direct_styles[i], label=solver)

                    plt.xlabel('Time $[s]$')
                    plt.ylabel('Norm of Residual')
                    plt.yscale('log')
                    plt.gca().set_ylim(top=10.0, bottom=1.0e-14)
                    if sort_system and preconditioner == "DIAG":
                        plt.legend()
                    plt.savefig(plot_dir+"{0}{1}_{2}_{3}_accuracy_over_time.pdf".format(case_name, "sorted" if sort_system else "unsorted", preconditioner, refinement))
                    plt.savefig(plot_dir+"{0}{1}_{2}_{3}_accuracy_over_time.svg".format(case_name, "sorted" if sort_system else "unsorted", preconditioner, refinement))
                    plt.close()
