import json
import subprocess as sp
from tempfile import NamedTemporaryFile
import shutil
import csv
import numpy as np


mesh_dir = "studies/matrix_solvers/meshes/"
report_dir = "studies/matrix_solvers/reports/"
result_dir = "studies/matrix_solvers/results/"
data_dir = "studies/matrix_solvers/data/"
iteration_dir = "studies/matrix_solvers/iterations/"
json_ext = ".json"
vtk_ext = ".vtk"
stl_ext = ".stl"


def get_matrix_condition(input_filename):
    # Gets the full method runtime, norm of final residual, and the total iterations for the given input

    # Run
    sp.run(["./machline.exe", input_filename])

    # Read in matrix
    A = np.genfromtxt("A_mat.txt")

    # Get SVD
    S = np.linalg.svd(A, compute_uv=False)
    S_max = S[0]
    S_min = S[-1]
    cond = np.sqrt(S_max/S_min)

    return S_max, S_min, cond


def write_input_file(input_filename, mesh_root_name, v_inf, M, solver, refinement, preconditioner, sort_system, vtk_mesh, iter_file='none', mirror_plane='none'):
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
            "mirror_about" : mirror_plane
        },
        "solver" : {
            "matrix_solver" : solver,
            "preconditioner" : preconditioner,
            "sort_system" : sort_system,
            "write_A_and_b" : True,
            "iterative_solver_output" : iter_file
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


def get_condition_numbers(mesh_root_name, v_inf, M, mirror_plane, vtk_mesh):
    # Determines all the condition numbers
    # Assumes the mesh has 'coarse', 'medium', and 'fine' refinements available

    # Options to iterate through
    #solver_options = ["LU", "BJAC", "BSSOR", "QRUP", "FQRUP", "GMRES"]
    #solver_options = ["BJAC", "BSSOR"]
    #solver_options = ["BSSOR"]
    solver = "GMRES"
    #refinement_options = ["fine"]
    refinement_options = ["coarse", "medium", "fine"]
    preconditioner_options = ["DIAG", "none"]
    sort_system_options = [False, True]
    #sort_system_options = [True]
    N_avg = 1

    # We'll just overwrite the input every time
    input_filename = "studies/matrix_solvers/input.json"

    # Iterate
    data = []
    for sort_system in sort_system_options:
        for preconditioner in preconditioner_options:
            for refinement in refinement_options:

                # Create input file
                write_input_file(input_filename, mesh_root_name, v_inf, M, solver, refinement, preconditioner, sort_system, vtk_mesh, mirror_plane=mirror_plane)

                # Get run times
                result = get_matrix_condition(input_filename)

                # Store
                data.append([sort_system, preconditioner, refinement]+list(result))

    # Write data
    with open("studies/matrix_solvers/"+mesh_root_name+"condition_data.csv", 'w') as data_handle:

        # header
        print("Sort System,Preconditioner,Refinement,S_max,S_min,Condition Number", file=data_handle)

        for data_item in data:
            print(",".join([str(x) for x in data_item]), file=data_handle)


if __name__=="__main__":

    # Run for cone
    case_root_name = "cone_10_deg_" # [coarse = 321, medium = 1241, fine = 4881]
    get_condition_numbers(case_root_name, [-1.0, 0.0, 0.0], 1.5, 'xy', True)

    # Run for diamond wing [coarse = 800, medium = 2152, fine = 7535]
    case_root_name = "diamond_5_deg_full_"
    get_condition_numbers(case_root_name, [1.0, 0.0, 0.0], 2.0, 'None', False)

    # Run for full configuration [coarse = 600, medium = 2400, fine = 11000]
    case_root_name = "full_config_"
    get_condition_numbers(case_root_name, [1.0, 0.0, 0.0], 2.0, 'xz', False)
