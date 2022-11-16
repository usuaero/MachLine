import os
import json
import numpy as np
import subprocess as sp
from tempfile import NamedTemporaryFile
import shutil
import csv


mesh_dir = "studies/matrix_solvers/meshes/"
report_dir = "studies/matrix_solvers/reports/"
result_dir = "studies/matrix_solvers/results/"
data_dir = "studies/matrix_solvers/data/"
iteration_dir = "studies/matrix_solvers/iterations/"
json_ext = ".json"
vtk_ext = ".vtk"
stl_ext = ".stl"


def get_full_method_results(input_filename):
    # Gets the full method runtime, norm of final residual, and the total iterations for the given input

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
    runtime = report["total_runtime"]
    if solver_stat == 0:
        sort_time = report["solver_results"]["timing"]["system_sorting"]
        prec_time = report["solver_results"]["timing"]["preconditioner"]
        solver_time = report["solver_results"]["timing"]["matrix_solver"]
        res_norm = report["solver_results"]["residual"]["norm"]
        iterations = report["solver_results"].get("iterations", "N/A")
        return runtime, sort_time, prec_time, solver_time, res_norm, iterations
    else:
        return runtime, solver_stat


def get_solver_runtime(solver):
    # Gets the runtime for the given solver
    
    # Run
    result = sp.run(["./solver_timer.exe", solver, 'A_mat.txt', 'b_vec.txt'], capture_output=True, text=True)
    if "Solution" in result.stdout:
        print(result.stdout)
        runtime = float(result.stdout.split()[-2])
        return runtime
    else:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError("Solver timer failed.")


def write_input_file(input_filename, mesh_root_name, v_inf, M, solver, refinement, preconditioner, sort_system, write_system, vtk_mesh, iter_file='none', mirror_plane='none'):
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
            "write_A_and_b" : write_system,
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


def run_paces(mesh_root_name, v_inf, M, mirror_plane, vtk_mesh):
    # Runs the MachLine matrix solvers through their paces 
    # Assumes the mesh has 'coarse', 'medium', and 'fine' refinements available

    # Options to iterate through
    solver_options = ["LU", "BJAC", "BSSOR", "QRUP", "FQRUP", "GMRES"]
    #solver_options = ["BJAC", "BSSOR"]
    #solver_options = ["BSSOR"]
    #refinement_options = ["fine"]
    refinement_options = ["coarse", "medium", "fine"]
    preconditioner_options = ["DIAG", "none"]
    sort_system_options = [False, True]
    #sort_system_options = [True]
    N_avg = 5

    # We'll just overwrite the input every time
    input_filename = "studies/matrix_solvers/input.json"

    # Iterate
    for sort_system in sort_system_options:
        for solver in solver_options:
            for preconditioner in preconditioner_options:
                for refinement in refinement_options:

                    # If this is an iterative solver, run one time writing out the residual history
                    if solver in ["GMRES", "BJAC", "BSSOR"]:
                        iter_file = iteration_dir + "{0}{1}_{2}_{3}_{4}_prec_history.csv".format(mesh_root_name, solver, refinement, preconditioner, "sorted" if sort_system else "unsorted")
                        write_input_file(input_filename, mesh_root_name, v_inf, M, solver, refinement, preconditioner, sort_system, False, vtk_mesh, iter_file=iter_file, mirror_plane=mirror_plane)
                        sp.run(["./machline.exe", input_filename])

                    for i in range(N_avg):

                        # Create input file
                        write_input_file(input_filename, mesh_root_name, v_inf, M, solver, refinement, preconditioner, sort_system, False, vtk_mesh, mirror_plane=mirror_plane)

                        # Get run times
                        result = get_full_method_results(input_filename)
                        if len(result) == 6:
                            runtime_m, sort_time, prec_time, solver_time, res_norm, iterations = result

                            # Write results
                            update_csv(mesh_root_name, solver, refinement, preconditioner, sort_system, i, runtime_m, solver_time+sort_time+prec_time, res_norm, iterations)

                        # Execution failed
                        else:

                            # Write results
                            update_csv(mesh_root_name, solver, refinement, preconditioner, sort_system, i, 'N/A', 'N/A', 'N/A', 'N/A')


def update_csv(mesh_root_name, solver, refinement, preconditioner, sort_system, i, runtime_m, runtime_s, res_norm, iterations):
    # Updates the csv data file

    # Initialize temporary file
    filename = "studies/matrix_solvers/data/"+mesh_root_name+"solver_test_data.csv"
    tempfile = NamedTemporaryFile(mode='w', delete=False)

    # Get column names
    fields = ["Solver", "Mesh Refinement", "Preconditioner", "Sorted", "Trial", "Method Run Time", "Solver Run Time", "Norm of Final Residual", "Iterations"]

    # Open file
    with open(filename, 'r') as csvfile, tempfile:

        # Get reader and writer
        reader = csv.DictReader(csvfile, fieldnames=fields)
        writer = csv.DictWriter(tempfile, fieldnames=fields)
        
        # Loop through rows to find the one we need
        for row in reader:

            # Update this row
            if row["Solver"] == solver and row["Mesh Refinement"] == refinement and row["Preconditioner"] == preconditioner and row["Sorted"].title() == str(sort_system) and int(row["Trial"]) == i:

                # Update row
                row["Method Run Time"] = runtime_m
                row["Solver Run Time"] = runtime_s
                row["Norm of Final Residual"] = res_norm
                row["Iterations"] = iterations

            # Write
            writer.writerow(row)

    # Move file
    shutil.move(tempfile.name, filename)


if __name__=="__main__":

    # Perform serial compilation of MachLine
    print("Compiling MachLine...")
    result = sp.run(["make", "serial"], capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError("MachLine compilation failed.")

    # Run for cone
    case_root_name = "cone_10_deg_" # [coarse = 321, medium = 1241, fine = 4881]
    run_paces(case_root_name, [-1.0, 0.0, 0.0], 1.5, 'xy', True)

    # Run for diamond wing [coarse = 800, medium = 2152, fine = 7535]
    case_root_name = "diamond_5_deg_full_"
    run_paces(case_root_name, [1.0, 0.0, 0.0], 2.0, 'None', False)

    # Run for full configuration [coarse = 600, medium = 2400, fine = 11000]
    case_root_name = "full_config_"
    run_paces(case_root_name, [1.0, 0.0, 0.0], 2.0, 'xz', False)
