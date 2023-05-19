import json
import subprocess as sp
import shutil
import numpy as np


mesh_dir = "studies/matrix_solvers/meshes/"
report_dir = "studies/matrix_solvers/reports/"
result_dir = "studies/matrix_solvers/results/"
data_dir = "studies/matrix_solvers/data/"
iteration_dir = "studies/matrix_solvers/iterations/"
json_ext = ".json"
vtk_ext = ".vtk"
stl_ext = ".stl"


def get_system_lower_bandwidth(input_filename):
    # Gets thelower bandwidth of the system for the given input file

    # Run
    sp.run(["./machline.exe", input_filename])

    # Get path to report file
    with open(input_filename, 'r') as input_handle:
        input_dict = json.load(input_handle)
        report_file = input_dict["output"]["report_file"]

    # Read in report
    with open(report_file, 'r') as report_handle:
        report = json.load(report_handle)

    # Check status
    solver_stat = report["solver_results"]["solver_status_code"]
    
    # Get results
    if solver_stat == 0:
        B_l = report["solver_results"]["system_lower_bandwidth"]
        N = report["solver_results"]["system_dimension"]
        return B_l, N
    else:
        raise RuntimeError("Solver failed.")


def write_input_file(input_filename, mesh_root_name, v_inf, M, order, refinement, preconditioner, sort_system, vtk_mesh, mirror_plane='none'):
    # Writes an input file

    # Assemble
    mesh_name = mesh_root_name + refinement
    if vtk_mesh:
        mesh_file = mesh_dir + mesh_name + vtk_ext
    else:
        mesh_file = mesh_dir + mesh_name + stl_ext
    input_dict ={
        "flow" : {
            "freestream_velocity" : list(v_inf),
            "freestream_mach_number" : M
        },
        "geometry" : {
            "file" : mesh_file,
            "spanwise_axis" : "+y",
            "mirror_about" : mirror_plane,
            "singularity_order" : order
        },
        "solver" : {
            "matrix_solver" : "FQRUP",
            "preconditioner" : preconditioner,
            "sort_system" : sort_system
        },
        "post_processing" : {
        },
        "output" : {
            "body_file" : result_dir + mesh_name + vtk_ext,
            "report_file" : report_dir + mesh_name + json_ext
        }
    }

    # Write
    with open(input_filename, 'w') as input_handle:
        json.dump(input_dict, input_handle, indent=4)


def analyze_lower_bandwidth(mesh_root_name, v_inf, M, mirror_plane, vtk_mesh):
    # Determines lower bandwidths for the different meshes using both lower- and higher-order distributions
    # Assumes the mesh has 'coarse', 'medium', and 'fine' refinements available

    # Options to iterate through
    refinement_options = ["coarse", "medium", "fine"]
    singularity_options = ["lower", "higher"]

    # We'll just overwrite the input every time
    input_filename = "studies/matrix_solvers/input.json"

    # Iterate
    Ns = np.zeros((3,2))
    B_ls = np.zeros((3,2))
    for j, order in enumerate(singularity_options):
        for i, refinement in enumerate(refinement_options):

            # Create input file
            write_input_file(input_filename, mesh_root_name, v_inf, M, order, refinement, "DIAG", True, vtk_mesh, mirror_plane=mirror_plane)

            # Get lower bandwidtn and system dimension
            B_ls[i,j], Ns[i,j] = get_system_lower_bandwidth(input_filename)

    print()
    print("         B_l")
    print(" N       Higher-Order   Lower-Order")
    print("------------------------------------------")
    for i in range(3):
        print(Ns[i,0], B_ls[i,1], B_ls[i,0])


if __name__=="__main__":

    # Run for cone
    case_root_name = "cone_10_deg_" # [coarse = 321, medium = 1241, fine = 4881]
    print("CONE")
    analyze_lower_bandwidth(case_root_name, [-1.0, 0.0, 0.0], 1.5, 'xy', True)

    # Run for diamond wing [coarse = 800, medium = 2152, fine = 7535]
    case_root_name = "diamond_5_deg_full_"
    print("DIAMOND WING")
    analyze_lower_bandwidth(case_root_name, [1.0, 0.0, 0.0], 2.0, 'None', False)

    # Run for full configuration [coarse = 600, medium = 2400, fine = 11000]
    case_root_name = "full_config_"
    print("FULL CONFIGURATION")
    analyze_lower_bandwidth(case_root_name, [1.0, 0.0, 0.0], 2.0, 'xz', False)
