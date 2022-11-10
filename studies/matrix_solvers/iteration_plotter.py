import numpy as np
import matplotlib.pyplot as plt


iteration_dir = "studies/matrix_solvers/iterations/"

# Options to iterate through
solver_options = ["BJAC", "BSSOR", "GMRES"]
refinement_options = ["coarse", "medium", "fine"]
preconditioner_options = ["DIAG", "none"]
sort_system_options = [False, True]


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

    # Test
    for solver in solver_options:
        for refinement in refinement_options:
            for preconditioner in preconditioner_options:
                for sort_system in sort_system_options:
                    res = get_residual_history("full_config_", solver, refinement, preconditioner, sort_system)
                    print(res)