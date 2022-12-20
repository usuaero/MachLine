import numpy as np
import matplotlib.pyplot as plt


iteration_dir = "studies/matrix_solvers/iterations/"
plot_dir = "studies/matrix_solvers/plots/"


def get_residual_history(mesh_root_name, solver, refinement, preconditioner, sort_system):
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


if __name__=="__main__":

    # Options
    case_names = ["cone_10_deg_", "diamond_5_deg_full_", "full_config_"]
    solver_options = ["BJAC", "BSSOR", "GMRES"]
    l_styles = ['-', '--', '-.']
    refinement_options = ["coarse", "medium", "fine"]
    colors = ['black', 'dimgray', 'silver']
    preconditioner_options = ["DIAG", "none"]
    sort_system_options = [False, True]

    # Run through
    for case_name in case_names:
        for sort_system in sort_system_options:
            for k, preconditioner in enumerate(preconditioner_options):
                plt.figure()
                for i, solver in enumerate(solver_options):
                    for j, refinement in enumerate(refinement_options):
                        res = get_residual_history(case_name, solver, refinement, preconditioner, sort_system)
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
