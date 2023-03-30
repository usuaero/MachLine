import os
import sys
import json

import subprocess as sp

from copy import deepcopy


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

    # Write file
    new_filename = original_filename.replace(".json", "_QUAD_{0}-order_{1}.json".format(order, formulation))
    write_input_file(copy_dict, new_filename)

    return new_filename


def run_machline(input_filename, delete_input=True):
    """Runs MachLine with the given input and returns the report if MachLine generated one."""

    # Run
    sp.run(["./machline.exe", input_filename])

    # Get report
    with open(input_filename, 'r') as input_handle:
        input_dict = json.load(input_handle)
    report_file = input_dict["output"].get("report_file")
    if report_file is not None:
        with open(report_file) as report_handle:
            report = json.load(report_handle)
    else:
        report = None

    # Delete input
    if delete_input:
        os.remove(input_filename)

    return report


def run_quad(input_filename, delete_input=True):
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
            reports.append(run_machline(altered_input_filename, delete_input=delete_input))

    return reports