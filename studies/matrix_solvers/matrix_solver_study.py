import os
import json
import numpy as np
import subprocess as sp


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
        solver_time = report["solver_results"]["matrix_solver_time"]
        res_norm = report["solver_results"]["residual"]["norm"]
        iterations = report["solver_results"].get("iterations", "N/A")
        return runtime, solver_time, res_norm, iterations
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
    # Assumes the mesh has 'coase', 'medium', and 'fine' refinements available

    # Options to iterate through
    solver_options = ["QRUP", "BJAC", "BSOR", "LU", "FQRUP", "GMRES"]
    refinement_options = ["coarse", "medium", "fine"]
    preconditioner_options = ["DIAG", "none"]
    sort_system_options = [True, False]
    N_avg = 5

    # We'll just overwrite the input every time
    input_filename = "studies/matrix_solvers/input.json"

    # Start up output
    with open("studies/matrix_solvers/data/"+mesh_root_name+"solver_test_data.csv", 'w') as data_handle:

        # Write out header
        print("Solver,Mesh Refinement,Preconditioner,Sorted,Trial,Method Run Time,Solver Run Time,Norm of Final Residual,Iterations", file=data_handle)
        data_handle.flush()

        # Iterate
        for solver in solver_options:
            for refinement in refinement_options:
                for preconditioner in preconditioner_options:
                    for sort_system in sort_system_options:

                        # Run MachLine once to write out linear system (need to set fast solver here)
                        write_input_file(input_filename, mesh_root_name, v_inf, M, "GMRES", refinement, preconditioner, sort_system, True, vtk_mesh, mirror_plane=mirror_plane)
                        sp.run(["./machline.exe", input_filename])

                        # If this is an iterative solver, run one time writing out the residual history
                        if solver in ["GMRES", "BJAC", "BSOR"]:
                            iter_file = iteration_dir + "{0}{1}_{2}_{3}_{4}_prec_history.csv".format(mesh_root_name, solver, refinement, preconditioner, "sorted" if sort_system else "unsorted")
                            write_input_file(input_filename, mesh_root_name, v_inf, M, solver, refinement, preconditioner, sort_system, False, vtk_mesh, iter_file=iter_file, mirror_plane=mirror_plane)
                            sp.run(["./machline.exe", input_filename])

                        for i in range(N_avg):

                            # Create input file
                            write_input_file(input_filename, mesh_root_name, v_inf, M, solver, refinement, preconditioner, sort_system, False, vtk_mesh, mirror_plane=mirror_plane)

                            # Get run times
                            result = get_full_method_results(input_filename)
                            if len(result) == 4:
                                runtime_m, runtime_s, res_norm, iterations = result

                                # Write results
                                print("{0},{1},{2},{3},{4},{5:.12e},{6:.12e},{7:.12e},{8}".format(solver, refinement, preconditioner, sort_system, i, runtime_m, runtime_s, res_norm, iterations), file=data_handle)
                                data_handle.flush()

                            # Execution failed
                            else:

                                # Write results
                                print("{0},{1},{2},{3},{4},{5},{6},{7},{8}".format(solver, refinement, preconditioner, sort_system, i, 'N/A', 'N/A', 'N/A', 'N/A'), file=data_handle)
                                data_handle.flush()


if __name__=="__main__":

    ## Compile solver timer
    #print("Compiling timer...")
    #result = sp.run(["gfortran", "-O2", "-fbounds-check", "-fbacktrace", "-fdefault-real-8", "common/linalg.f95", "dev/time_matrix_solver.f95", "-o", "solver_timer.exe"], capture_output=True, text=True)
    #if result.returncode != 0:
    #    print(result.stdout)
    #    print(result.stderr)
    #    raise RuntimeError("Solver timer compilation failed.")

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
