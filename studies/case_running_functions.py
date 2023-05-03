import os
import json
import numpy as np

import subprocess as sp
import multiprocessing as mp

from copy import deepcopy


cases = ['ML', 'MH', 'SL', 'SH']
line_styles = ['k-', 'k--', 'k:', 'k-.']
quad_labels = ['_QUAD_lower-order_morino', '_QUAD_higher-order_morino', '_QUAD_lower-order_source-free', '_QUAD_higher-order_source-free', ]


def write_input_file(input_dict, input_filename):
    """Writes the given input dict to the given file location."""

    with open(input_filename, 'w') as input_handle:
        json.dump(input_dict, input_handle, indent=4)


def write_altered_input_file(original_filename, input_dict, order, formulation):
    """Updates the input dict and writes a new file. For use with run_quad()."""

    # Create copy
    copy_dict = deepcopy(input_dict)
    
    # Update options
    copy_dict["solver"]["formulation"] = formulation
    copy_dict["geometry"]["singularity_order"] = order

    # Update outputs
    for key, value in input_dict["output"].items():
        if "vtk" in value:
            copy_dict["output"][key] = value.replace(".vtk", "_QUAD_{0}-order_{1}.vtk".format(order, formulation))
        elif "json" in value:
            copy_dict["output"][key] = value.replace(".json", "_QUAD_{0}-order_{1}.json".format(order, formulation))

    # Write file
    new_filename = original_filename.replace(".json", "_QUAD_{0}-order_{1}.json".format(order, formulation))
    write_input_file(copy_dict, new_filename)

    return new_filename


def run_machline(input_filename, delete_input=True, run=True):
    """Runs MachLine with the given input and returns the report if MachLine generated one."""

    # Run
    if run:
        sp.run(["./machline.exe", input_filename])

    # Get report
    with open(input_filename, 'r') as input_handle:
        input_dict = json.load(input_handle)
    report_file = input_dict["output"].get("report_file")
    if report_file is not None:
        try:
            with open(report_file) as report_handle:
                report = json.load(report_handle)
        except:
            report = None
    else:
        report = None

    # Delete input
    if delete_input:
        os.remove(input_filename)

    return report


def _run_order_and_formulation(*args):
    # Runs the given order and formulation

    input_filename = args[0]
    input_dict = args[1]
    order = args[2]
    formulation = args[3]
    delete_input = args[4]
    run = args[5]

    # Write out new input
    altered_input_filename = write_altered_input_file(input_filename, input_dict, order, formulation)

    # Run MachLine
    return run_machline(altered_input_filename, delete_input=delete_input, run=run)


def run_quad_parallel(input_filename, delete_input=True, run=True):
    """Runs the given filename with the four option combinations in parallel

    The reports will be returned in the order:
        1-Morino, lower-order
        2-Morino, higher-order
        3-Source-free, lower-order
        4-Source-free, higher-order
    """

    # Load input
    with open(input_filename, 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Options
    formulations = ["morino", "source-free"]
    orders = ["lower", "higher"]

    # Get inputs and run
    arg_list = []
    for formulation in formulations:
        for order in orders:
            arg_list.append((input_filename, input_dict, order, formulation, delete_input, run))

    with mp.Pool() as pool:
        reports = pool.map(_run_order_and_formulation, arg_list)

    return reports


def run_quad(input_filename, delete_input=True, run=True):
    """Runs the given filename with the four option combinations

    The reports will be returned in the order:
        1-Morino, lower-order
        2-Morino, higher-order
        3-Source-free, lower-order
        4-Source-free, higher-order
    """

    # Load input
    with open(input_filename, 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Options
    formulations = ["morino", "source-free"]
    orders = ["lower", "higher"]

    # Get inputs and run
    reports = []
    for formulation in formulations:
        for order in orders:

            # Write out new input
            altered_input_filename = write_altered_input_file(input_filename, input_dict, order, formulation)

            # Run MachLine
            reports.append(run_machline(altered_input_filename, delete_input=delete_input, run=run))

    return reports


def get_order_of_convergence(grid_parameter, result, truth_from_results=True):
    """Returns the order of convergence from the provided results.
    
    Parameters
    ----------
    grid_parameter : list or ndarray
        List of grid parameters, such as number of vertices, length scales, against which the convergence will be measured.

    result : list or ndarray
        The resulting value of which the convergence is to be determined. Should be the same length as grid_parameter.

    truth_from_results : bool, optional
        Whether the truth result should be determined from the results using Richarson extrapolation. Defaults to True.
        If set to False, the truth is assumed to be zero.
    """

    # Get the error
    if truth_from_results:
        err = np.abs((result[:-1] - result[-1])/result[-1])
        coefs = np.polyfit(np.log(grid_parameter[:-1]), np.log(err), deg=1)
    else:
        err = np.abs(result)
        coefs = np.polyfit(np.log(grid_parameter), np.log(err), deg=1)

    # Calculate order of convergence
    return np.abs(coefs[0].item())