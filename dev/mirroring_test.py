# Tests the matching achieved by mirroring

import json

import subprocess as sp
import numpy as np


def get_min_max_cp():
    # Gets the min and max pressure coefficients from the report file

    with open('dev/report.txt', 'r') as report_file:

        lines = report_file.readlines()
        C_p_max = float(lines[3].split()[-1])
        C_p_min = float(lines[4].split()[-1])

    return C_p_max, C_p_min


if __name__=="__main__":

    # Initialize
    full_input_file = 'dev/full_wing_input.json'
    half_input_file = 'dev/half_wing_input.json'
    exe = './machline.exe'

    # Options
    flow_conditions = ['sym', 'asym']
    waked_options = [False, True]

    # Begin table
    print()
    print("# of Decimal Places Matching in Pressure Coefficients")
    print("-----------------------------------------------------")
    print()
    print("".join([" "]*40)+"|  {0:<20}|  {1:<20}|  {2:<20}".format("w/ Full", "w/ Full, -V", "w/ Self, -V"))
    print("".join(["-"]*110))

    # Load inputs
    with open(full_input_file, 'r') as full_input_handle:
        full_input = json.load(full_input_handle)
    with open(half_input_file, 'r') as half_input_handle:
        half_input = json.load(half_input_handle)

    # Loop through flow conditions
    for flow_condition in flow_conditions:

        # Print out flow condition
        print("{0:>15} flow".format(flow_condition)+"".join([" "]*20)+"|"+"".join([" "]*22)+"|"+"".join([" "]*22)+"|")

        # Set flow velocity
        if flow_condition == 'sym':
            full_input["flow"]["freestream_velocity"] = [100.0, 0.0, 10.0]
            half_input["flow"]["freestream_velocity"] = [100.0, 0.0, 10.0]
        else:
            full_input["flow"]["freestream_velocity"] = [100.0, 30.0, 10.0]
            half_input["flow"]["freestream_velocity"] = [100.0, 30.0, 10.0]

        # Loop through wake options
        for waked in waked_options:

            # Has wake
            if waked:

                # Set wake
                full_input["geometry"]["wake_model"]["wake_shedding_angle"] = 90.0
                half_input["geometry"]["wake_model"]["wake_shedding_angle"] = 90.0

                # Print out info
                print("".join([" "]*20)+"              waked "+"|"+"".join([" "]*22)+"|"+"".join([" "]*22)+"|")

            # No wake
            else:

                # Set wake
                full_input["geometry"]["wake_model"]["wake_shedding_angle"] = 179.0
                half_input["geometry"]["wake_model"]["wake_shedding_angle"] = 179.0

                # Print out info
                print("".join([" "]*20)+"           wakeless "+"|"+"".join([" "]*22)+"|"+"".join([" "]*22)+"|")

            # Compare full mesh to half mesh with original velocity

            # Write out inputs
            with open(full_input_file, 'w') as full_input_handle:
                json.dump(full_input, full_input_handle, indent=4)
            with open(half_input_file, 'w') as half_input_handle:
                json.dump(half_input, half_input_handle, indent=4)

            # Run full mesh
            with sp.Popen([exe, full_input_file], stdin=sp.PIPE, stdout=sp.PIPE) as mftran_process:
                mftran_process.wait()

            # Get pressure coefs
            C_p_max_full, C_p_min_full = get_min_max_cp()

            # Run half mesh
            with sp.Popen([exe, half_input_file], stdin=sp.PIPE, stdout=sp.PIPE) as mftran_process:
                mftran_process.wait()

            # Get pressure coefs
            C_p_max_half, C_p_min_half = get_min_max_cp()

            # Calculate maximum difference
            max_diff = max(abs(C_p_max_full-C_p_max_half), abs(C_p_min_full-C_p_min_half))

            # Number of matching decimal places
            N_match_full_half = -round(np.log10(max_diff))

            # Compare full mesh to half mesh with flipped velocity

            # Set flow velocity
            if flow_condition == 'sym':
                full_input["flow"]["freestream_velocity"] = [100.0, 0.0, -10.0]
                half_input["flow"]["freestream_velocity"] = [100.0, 0.0, -10.0]
            else:
                full_input["flow"]["freestream_velocity"] = [100.0, -30.0, 10.0]
                half_input["flow"]["freestream_velocity"] = [100.0, -30.0, 10.0]

            # Write out inputs
            with open(full_input_file, 'w') as full_input_handle:
                json.dump(full_input, full_input_handle, indent=4)
            with open(half_input_file, 'w') as half_input_handle:
                json.dump(half_input, half_input_handle, indent=4)

            # Run full mesh
            with sp.Popen([exe, full_input_file], stdin=sp.PIPE, stdout=sp.PIPE) as mftran_process:
                mftran_process.wait()

            # Get pressure coefs
            C_p_max_full, C_p_min_full = get_min_max_cp()

            # Run half mesh
            with sp.Popen([exe, half_input_file], stdin=sp.PIPE, stdout=sp.PIPE) as mftran_process:
                mftran_process.wait()

            # Get pressure coefs
            C_p_max_half, C_p_min_half = get_min_max_cp()

            # Calculate maximum difference
            max_diff = max(abs(C_p_max_full-C_p_max_half), abs(C_p_min_full-C_p_min_half))

            # Number of matching decimal places
            N_match_full_half_flipped = -round(np.log10(max_diff))

            # Compare half mesh to self with flipped velocity

            # Set flow velocity
            if flow_condition == 'sym':
                half_input["flow"]["freestream_velocity"] = [100.0, 0.0, 10.0]
            else:
                half_input["flow"]["freestream_velocity"] = [100.0, 30.0, 10.0]

            # Write out input
            with open(half_input_file, 'w') as half_input_handle:
                json.dump(half_input, half_input_handle, indent=4)

            # Run half mesh
            with sp.Popen([exe, half_input_file], stdin=sp.PIPE, stdout=sp.PIPE) as mftran_process:
                mftran_process.wait()

            # Get pressure coefs
            C_p_max_orig, C_p_min_orig = get_min_max_cp()

            # Set flow velocity
            if flow_condition == 'sym':
                half_input["flow"]["freestream_velocity"] = [100.0, 0.0, -10.0]
            else:
                half_input["flow"]["freestream_velocity"] = [100.0, -30.0, 10.0]

            # Write out input
            with open(half_input_file, 'w') as half_input_handle:
                json.dump(half_input, half_input_handle, indent=4)

            # Run half mesh
            with sp.Popen([exe, half_input_file], stdin=sp.PIPE, stdout=sp.PIPE) as mftran_process:
                mftran_process.wait()

            # Get pressure coefs
            C_p_max_flip, C_p_min_flip = get_min_max_cp()

            # Calculate maximum difference
            max_diff = max(abs(C_p_max_orig-C_p_max_flip), abs(C_p_min_orig-C_p_min_flip))

            # Number of matching decimal places
            N_match_half_half_flipped = -round(np.log10(max_diff))

            print("".join([" "]*40)+"|  {0:<20}|  {1:<20}|  {2:<20}".format(N_match_full_half, N_match_full_half_flipped, N_match_half_half_flipped))
            print("                    "+"".join(["-"]*90))

        print("".join(["-"]*110))