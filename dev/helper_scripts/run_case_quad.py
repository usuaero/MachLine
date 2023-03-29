import os
import sys
import json

import subprocess as sp

from copy import deepcopy


def write_altered_input_file(original_filename, input_dict, order, formulation):
    # Updates the input dict and writes a new file

    # Create copy
    copy_dict = deepcopy(input_dict)
    
    # Update options
    copy_dict["solver"]["formulation"] = formulation
    copy_dict["geometry"]["singularity_order"] = order

    # Update outputs
    for key, value in input_dict["output"].items():
        if "vtk" in value:
            copy_dict["output"][key] = value.replace(".vtk", "_QUAD_{0}-order_{1}.vtk".format(order, formulation))

    # Get new filename
    new_filename = original_filename.replace(".json", "_QUAD_{0}-order_{1}.json".format(order, formulation))

    # Dump
    with open(new_filename, 'w') as new_handle:
        json.dump(copy_dict, new_handle, indent=4)

    return new_filename


def run_machline(input_filename, delete_input=True):
    # Runs MachLine with the given input

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
    # Runs the given filename with the four option combinations

    # Load input
    with open(input_filename, 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Options
    formulations = ["source-free", "morino"]
    orders = ["lower", "higher"]

    # Get inputs and run
    for formulation in formulations:
        for order in orders:

            # Write out new input
            altered_input_filename = write_altered_input_file(input_filename, input_dict, order, formulation)

            # Run MachLine
            run_machline(altered_input_filename, delete_input=False)


if __name__=="__main__":

    # Get input file
    input_filename = sys.argv[-1]

    # Run
    run_quad(input_filename)