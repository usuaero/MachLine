import os
import json
import shutil
import subprocess as sp


class MachLineError(Exception):
    pass


def run_machline(input_file, remove_input=False):
    # Runs MachLine and delivers the output; cleans up output files
    
    # Create results directory
    if not os.path.exists("test/results"):
        os.mkdir("test/results")

    # Run MachLine
    result = sp.run(["./machline.exe", input_file], capture_output=True, text=True)

    # Remove input
    if remove_input:
        os.remove(input_file)

    # Read in report
    if os.path.exists("test/results/report.txt"):
        with open("test/results/report.txt", 'r') as report_handle:
            report = report_handle.readlines()

        # Check if MachLine thinks it was successful
        success = result.stdout.split()[-1] == "successfully."

    else:

        # If no report file was generated, MachLine was not successful
        success = False

    # Clean up results directory
    shutil.rmtree("test/results")

    if not success:
        print(result.stdout)
        raise MachLineError

    else:
        C_p_max = float(report[3].split()[-1])
        C_p_min = float(report[4].split()[-1])
        Cx = float(report[5].split()[-1])
        Cy = float(report[6].split()[-1])
        Cz = float(report[7].split()[-1])
        return C_p_max, C_p_min, Cx, Cy, Cz


def test_half_wing_source_free():
    # Tests the half wing case with the source-free formulation returns the consistent result

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/half_wing_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.6128404365547362) < 1e-12)
    assert(abs(C_p_min - -0.7562621922470163) < 1e-12)
    assert(abs(Cx - -0.006929293635976442) < 1e-12)
    assert(abs(Cy - -0.00024116090946916353) < 1e-12)
    assert(abs(Cz - -0.22240033375591964) < 1e-12)


def test_sphere_morino():
    # Tests the sphere case with the Morino formulation returns the consistent result

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/sphere_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.991141990209599) < 1e-12)
    assert(abs(C_p_min - -1.2378311776000124) < 1e-12)
    assert(abs(Cx - 1.2420239657237597e-07) < 1e-12)
    assert(abs(Cy - -1.557722467282474e-05) < 1e-12)
    assert(abs(Cz - -2.0422563784538994e-06) < 1e-12)


def test_full_half_wing_compare_morino():
    # Tests the full wing and half wing with the Morino formulation return identical results

    # Set full wing input
    with open("test/full_wing_input.json", 'r') as full_wing_input_handle:
        input_dict = json.load(full_wing_input_handle)

    input_dict["flow"]["freestream_velocity"] = [100.0, 0.0, 10.0]
    input_dict["solver"]["formulation"] = 'morino'

    with open("test/altered_full_wing_input.json", 'w') as full_wing_input_handle:
        input_dict = json.dump(input_dict, full_wing_input_handle, indent=4)

    # Run full wing
    C_p_max_full, C_p_min_full, Cx_full, Cy_full, Cz_full = run_machline("test/altered_full_wing_input.json", remove_input=True)

    # Set half wing input
    with open("test/half_wing_input.json", 'r') as half_wing_input_handle:
        input_dict = json.load(half_wing_input_handle)

    input_dict["flow"]["freestream_velocity"] = [100.0, 0.0, 10.0]
    input_dict["solver"]["formulation"] = 'morino'

    with open("test/altered_half_wing_input.json", 'w') as half_wing_input_handle:
        input_dict = json.dump(input_dict, half_wing_input_handle, indent=4)

    C_p_max_half, C_p_min_half, Cx_half, Cy_half, Cz_half = run_machline("test/altered_half_wing_input.json", remove_input=True)

    print(abs(C_p_max_full-C_p_max_half))
    print(abs(C_p_min_full-C_p_min_half))

    assert(abs(C_p_max_full-C_p_max_half)<1e-4)
    assert(abs(C_p_min_full-C_p_min_half)<2e-4)