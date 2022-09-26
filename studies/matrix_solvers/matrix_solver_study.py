import os
import json
import numpy as np
import subprocess as sp


def get_full_method_results(input_filename):
    # Gets the full method runtime, norm of final residual, and the total iterations for the given input

    # Run
    sp.run(["./machline.exe", input_filename])

    # Get results
    with open("studies/matrix_solvers/results/report.json", 'r') as report_handle:
        report = json.load(report_handle)
    runtime = report["total_runtime"]
    res_norm = report["solver_results"]["residual"]["norm"]
    iterations = report["solver_results"].get("iterations", "N/A")
    return runtime, res_norm, iterations


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

def write_input_file(input_filename, mesh_root_name, v_inf, M, solver, refinement, preconditioner, sort_system, write_system, iter_file='none', mirror_plane='none'):
    # Writes an input file

    # Assemble
    mesh_name = "studies/matrix_solvers/meshes/"+mesh_root_name+refinement+".vtk"
    input_dict ={
        "flow" : {
            "freestream_velocity" : list(v_inf),
            "freestream_mach_number" : M
        },
        "geometry" : {
            "file" : mesh_name,
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
            "body_file" : mesh_name.replace("meshes", "results"),
            "report_file" : "studies/matrix_solvers/results/report.json"
        }
    }

    # Write
    with open(input_filename, 'w') as input_handle:
        json.dump(input_dict, input_handle, indent=4)


def run_paces(mesh_root_name, v_inf, M, mirror_plane):
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
    with open("studies/matrix_solvers/results/"+mesh_root_name+"solver_test_data.csv", 'w') as data_handle:

        # Write out header
        print("Solver,Mesh Refinement,Preconditioner,Sorted,Trial,Method Run Time,Solver Run Time,Norm of Final Residual,Iterations", file=data_handle)
        data_handle.flush()

        # Iterate
        for solver in solver_options:
            for refinement in refinement_options:
                for preconditioner in preconditioner_options:
                    for sort_system in sort_system_options:

                        # Run MachLine once to write out linear system (need to set fast solver here)
                        write_input_file(input_filename, mesh_root_name, v_inf, M, "GMRES", refinement, preconditioner, sort_system, True, mirror_plane=mirror_plane)
                        sp.run(["./machline.exe", input_filename])

                        # If this is an iterative solver, run one time writing out the residual history
                        if solver in ["GMRES", "BJAC", "BSOR"]:
                            iter_file = "studies/matrix_solvers/results/{0}_{1}_{2}_prec_history.csv".format(solver, refinement, preconditioner)
                            write_input_file(input_filename, mesh_root_name, v_inf, M, solver, refinement, preconditioner, sort_system, False, iter_file=iter_file, mirror_plane=mirror_plane)
                            sp.run(["./machline.exe", input_filename])

                        for i in range(N_avg):

                            # Create input file
                            write_input_file(input_filename, mesh_root_name, v_inf, M, solver, refinement, preconditioner, sort_system, False, mirror_plane=mirror_plane)

                            # Get run time for all of MachLine
                            runtime_m, res_norm, iterations = get_full_method_results(input_filename)

                            # Get run time for matrix solver by itself
                            runtime_s = get_solver_runtime(solver)

                            # Write results
                            print("{0},{1},{2},{3},{4},{5},{6},{7},{8}".format(solver, refinement, preconditioner, sort_system, i, runtime_m, runtime_s, res_norm, iterations), file=data_handle)
                            data_handle.flush()




if __name__=="__main__":

    # Compile solver timer
    print("Compiling timer...")
    result = sp.run(["gfortran", "-O2", "-fbounds-check", "-fbacktrace", "-fdefault-real-8", "common/linalg.f95", "dev/time_matrix_solver.f95", "-o", "solver_timer.exe"], capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError("Solver timer compilation failed.")

    # Perform serial compilation of MachLine
    print("Compiling MachLine...")
    result = sp.run(["make", "serial"], capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError("MachLine compilation failed.")

    # Run for cone
    case_root_name = "cone_10_deg_"
    run_paces(case_root_name, [-1.0, 0.0, 0.0], 1.5, 'xy')