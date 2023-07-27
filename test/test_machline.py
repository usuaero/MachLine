# Regression tests for MachLine
# Please number each test using the format already in place
# KEEP THE TESTS IN ORDER!


import os
import json
import shutil
import subprocess as sp
import numpy as np

UPDATE_OFFBODY_RESULTS = False


class MachLineError(Exception):
    pass


def run_machline(input_file, remove_input=False, remove_results=True):
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
    report_file = "test/results/report.json"
    if os.path.exists(report_file):
        with open(report_file, 'r') as report_handle:
            report = json.load(report_handle)

        # Check if MachLine thinks it was successful
        success = "MachLine exited successfully." in result.stdout

    else:

        # If no report file was generated, MachLine was not successful
        success = False

    # Clean up results directory
    if remove_results:
        shutil.rmtree("test/results")

        # Other files to clean up
        if os.path.exists("iterative_solver_prog.csv"):
            os.remove("iterative_solver_prog.csv")

    if not success:
        print(result.stdout)
        print(result.stderr)
        raise MachLineError

    else:
        try:
            C_p_max = float(report['pressure_calculations']['incompressible_rule']['max'])
            C_p_min = float(report['pressure_calculations']['incompressible_rule']['min'])
        except KeyError:
            C_p_max = float(report['pressure_calculations']['isentropic_rule']['max'])
            C_p_min = float(report['pressure_calculations']['isentropic_rule']['min'])
        Cx = float(report['total_forces']['Cx'])
        Cy = float(report['total_forces']['Cy'])
        Cz = float(report['total_forces']['Cz'])
        return C_p_max, C_p_min, Cx, Cy, Cz


def test_01_half_wing_lower_source_free_asym_inc_flow():
    # Tests the half wing case with the source-free formulation returns the consistent result

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/half_wing_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.749427668724814) < 1e-10)
    assert(abs(C_p_min - -1.28122767724128) < 1e-9)
    assert(abs(Cx - -0.394147360647448) < 1e-9)
    assert(abs(Cy - -0.0469272618384146) < 1e-9)
    assert(abs(Cz - 20.6524631977163) < 1e-12)


def test_02_half_wing_higher_source_free_asym_inc_flow():
    # Tests the half wing case with the source-free formulation returns the consistent result

    # Load original input
    with open("test/input_files/half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Alter input
    input_dict["geometry"]["singularity_order"] = "higher"

    # Write altered input
    altered_input_file = "test/input_files/altered_half_wing_input.json"
    with open(altered_input_file, 'w') as altered_input_handle:
        json.dump(input_dict, altered_input_handle, indent=4)

    # Run MachLine
    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/altered_half_wing_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.749393430277945) < 1e-12)
    assert(abs(C_p_min - -1.1156939646846) < 1e-12)
    assert(abs(Cx - -0.3769033917615) < 1e-9)
    assert(abs(Cy - -0.0448030638527709) < 1e-9)
    assert(abs(Cz - 20.1961409335859) < 1e-12)


def test_03_half_wing_lower_morino_asym_inc_flow():
    # Tests the half wing case with the source-free formulation returns the consistent result

    # Load original input
    with open("test/input_files/half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Alter input
    input_dict["solver"]["formulation"] = "morino"

    # Write altered input
    altered_input_file = "test/input_files/altered_half_wing_input.json"
    with open(altered_input_file, 'w') as altered_input_handle:
        json.dump(input_dict, altered_input_handle, indent=4)

    # Run MachLine
    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/altered_half_wing_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.749426106253892) < 1e-10)
    assert(abs(C_p_min - -1.28120652757854) < 1e-9)
    assert(abs(Cx - -0.394142487065741) < 1e-9)
    assert(abs(Cy - -0.0469266124998036) < 1e-9)
    assert(abs(Cz - 20.6521815154214) < 1e-12)


def test_04_half_wing_higher_morino_asym_inc_flow():
    # Tests the half wing case with the source-free formulation returns the consistent result

    # Load original input
    with open("test/input_files/half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Alter input
    input_dict["solver"]["formulation"] = "morino"
    input_dict["geometry"]["singularity_order"] = "higher"

    # Write altered input
    altered_input_file = "test/input_files/altered_half_wing_input.json"
    with open(altered_input_file, 'w') as altered_input_handle:
        json.dump(input_dict, altered_input_handle, indent=4)

    # Run MachLine
    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/altered_half_wing_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.750968663253082) < 1e-10)
    assert(abs(C_p_min - -1.10238783762737) < 1e-10)
    assert(abs(Cx - -0.379595457329634) < 1e-9)
    assert(abs(Cy - -0.0459678458682901) < 1e-9)
    assert(abs(Cz - 20.2266509889293) < 1e-12)


def test_05_half_wing_lower_morino_zero_aoa_zero_beta_inc():
    # Tests the half wing case at zero angle of attack with the Morino formulation

    # Load original input
    with open("test/input_files/half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Alter input
    input_dict["flow"]["freestream_velocity"] = [100.0, 0.0, 0.0]
    input_dict["solver"]["formulation"] = "morino"

    # Write altered input
    altered_input_file = "test/input_files/altered_half_wing_input.json"
    with open(altered_input_file, 'w') as altered_input_handle:
        json.dump(input_dict, altered_input_handle, indent=4)

    # Run MachLine
    C_p_max, C_p_min, Cx, Cy, Cz = run_machline(altered_input_file, remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.221628441564136) < 1e-12)
    assert(abs(C_p_min - -0.427625486100214) < 1e-12)
    assert(abs(Cx - 0.301682907690379) < 1e-12)
    assert(abs(Cy) < 1e-12)
    assert(abs(Cz - 3.22319026937329e-12) < 1e-12)


def test_06_half_wing_higher_morino_zero_aoa_zero_beta_inc():
    # Tests the half wing case at zero angle of attack with the Morino formulation

    # Load original input
    with open("test/input_files/half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Alter input
    input_dict["flow"]["freestream_velocity"] = [100.0, 0.0, 0.0]
    input_dict["solver"]["formulation"] = "morino"
    input_dict["geometry"]["singularity_order"] = "higher"

    # Write altered input
    altered_input_file = "test/input_files/altered_half_wing_input.json"
    with open(altered_input_file, 'w') as altered_input_handle:
        json.dump(input_dict, altered_input_handle, indent=4)

    # Run MachLine
    C_p_max, C_p_min, Cx, Cy, Cz = run_machline(altered_input_file, remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.217923409759566) < 1e-12)
    assert(abs(C_p_min - -0.432426859568809) < 1e-12)
    assert(abs(Cx - 0.308691949491238) < 1e-12)
    assert(abs(Cy) < 1e-12)
    assert(abs(Cz) < 1e-11)


def test_07_sphere_morino_inc_iterative():
    # Tests the sphere case with the Morino formulation using an iterative matrix solver

    # Load original input
    with open("test/input_files/sphere_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Alter input
    input_dict["solver"]["matrix_solver"] = "BJAC"
    input_dict["solver"]["relaxation"] = 0.9

    # Write altered input
    altered_input_file = "test/input_files/altered_sphere_input.json"
    with open(altered_input_file, 'w') as altered_input_handle:
        json.dump(input_dict, altered_input_handle, indent=4)

    # Run MachLine
    C_p_max, C_p_min, Cx, Cy, Cz = run_machline(altered_input_file, remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.99114275830853) < 1e-12)
    assert(abs(C_p_min - -1.23782791839695) < 1e-12)
    assert(abs(Cx) < 2e-5)
    assert(abs(Cy) < 2e-5)
    assert(abs(Cz) < 2e-5)


def test_08_sphere_morino_inc():
    # Tests the sphere case with the Morino formulation returns the consistent result

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/sphere_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.991142758308531) < 1e-12)
    assert(abs(C_p_min + 1.23782791839695) < 1e-12)
    assert(abs(Cx) < 1e-4)
    assert(abs(Cy) < 1e-4)
    assert(abs(Cz) < 1e-4)


def test_09_full_half_wing_compare_lower_morino_zero_beta_inc():
    # Tests the full wing and half wing with the Morino formulation return identical results

    # Set full wing input
    with open("test/input_files/full_wing_input.json", 'r') as full_wing_input_handle:
        input_dict = json.load(full_wing_input_handle)

    input_dict["flow"]["freestream_velocity"] = [100.0, 0.0, 10.0]
    input_dict["solver"]["formulation"] = 'morino'

    with open("test/input_files/altered_full_wing_input.json", 'w') as full_wing_input_handle:
        input_dict = json.dump(input_dict, full_wing_input_handle, indent=4)

    # Run full wing
    C_p_max_full, C_p_min_full, Cx_full, Cy_full, Cz_full = run_machline("test/input_files/altered_full_wing_input.json", remove_input=True)

    # Set half wing input
    with open("test/input_files/half_wing_input.json", 'r') as half_wing_input_handle:
        input_dict = json.load(half_wing_input_handle)

    input_dict["flow"]["freestream_velocity"] = [100.0, 0.0, 10.0]
    input_dict["solver"]["formulation"] = 'morino'

    with open("test/input_files/altered_half_wing_input.json", 'w') as half_wing_input_handle:
        input_dict = json.dump(input_dict, half_wing_input_handle, indent=4)

    C_p_max_half, C_p_min_half, Cx_half, Cy_half, Cz_half = run_machline("test/input_files/altered_half_wing_input.json", remove_input=True)

    print(abs(C_p_max_full-C_p_max_half))
    print(abs(C_p_min_full-C_p_min_half))
    print(abs(Cx_full-Cx_half))
    print(abs(Cy_full-Cy_half))
    print(abs(Cz_full-Cz_half))

    assert(abs(C_p_max_full-C_p_max_half)<1e-3)
    assert(abs(C_p_min_full-C_p_min_half)<2e-3)
    assert(abs(Cx_full-Cx_half)<1e-3)
    assert(abs(Cy_full-Cy_half)<1e-3)
    assert(abs(Cz_full-Cz_half)<2e-2)


def test_10_full_half_wing_compare_higher_morino_zero_beta_inc():
    # Tests the full wing and half wing with the Morino formulation return identical results

    # Set full wing input
    with open("test/input_files/full_wing_input.json", 'r') as full_wing_input_handle:
        input_dict = json.load(full_wing_input_handle)

    input_dict["flow"]["freestream_velocity"] = [100.0, 0.0, 10.0]
    input_dict["solver"]["formulation"] = 'morino'
    input_dict["geometry"]["singularity_order"] = "higher"

    with open("test/input_files/altered_full_wing_input.json", 'w') as full_wing_input_handle:
        input_dict = json.dump(input_dict, full_wing_input_handle, indent=4)

    # Run full wing
    C_p_max_full, C_p_min_full, Cx_full, Cy_full, Cz_full = run_machline("test/input_files/altered_full_wing_input.json", remove_input=True)

    # Set half wing input
    with open("test/input_files/half_wing_input.json", 'r') as half_wing_input_handle:
        input_dict = json.load(half_wing_input_handle)

    input_dict["flow"]["freestream_velocity"] = [100.0, 0.0, 10.0]
    input_dict["solver"]["formulation"] = 'morino'
    input_dict["geometry"]["singularity_order"] = "higher"

    with open("test/input_files/altered_half_wing_input.json", 'w') as half_wing_input_handle:
        input_dict = json.dump(input_dict, half_wing_input_handle, indent=4)

    C_p_max_half, C_p_min_half, Cx_half, Cy_half, Cz_half = run_machline("test/input_files/altered_half_wing_input.json", remove_input=True)

    print(abs(C_p_max_full-C_p_max_half))
    print(abs(C_p_min_full-C_p_min_half))
    print(abs(Cx_full-Cx_half))
    print(abs(Cy_full-Cy_half))
    print(abs(Cz_full-Cz_half))

    assert(abs(C_p_max_full-C_p_max_half)<1e-3)
    assert(abs(C_p_min_full-C_p_min_half)<2e-3)
    assert(abs(Cx_full-Cx_half)<1e-3)
    assert(abs(Cy_full-Cy_half)<1e-3)
    assert(abs(Cz_full-Cz_half)<2e-2)


def test_11_subsonic_comp_pressure_corrections():

    # Tests the half wing case with the source-free formulation returns the consistent result

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/compressible_half_wing_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.817385304260662) < 1e-9)
    assert(abs(C_p_min - -1.09394053248534) < 1e-9)
    assert(abs(Cx - -0.508981218975681) < 1e-9)
    assert(abs(Cy - 0.0102074574010827) < 1e-9)
    assert(abs(Cz - 26.2860762052636) < 1e-8)


def test_12_supersonic_half_wing_morino_zero_aoa_zero_beta():

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/supersonic_half_wing_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.121697468024553) < 1e-12)
    assert(abs(C_p_min - -0.116094421618336) < 1e-12)
    assert(abs(Cx - 0.1424008894449) < 1e-12)
    assert(abs(Cy) < 1e-12)
    assert(abs(Cz) < 1e-12)


def test_13_supersonic_half_wing_lower_source_free_zero_aoa_zero_beta():

    with open("test/input_files/supersonic_half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    input_dict["solver"]["formulation"] = 'source-free'

    with open("test/input_files/altered_supersonic_input.json", 'w') as input_handle:
        input_dict = json.dump(input_dict, input_handle, indent=4)

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/altered_supersonic_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.121697468024034) < 1e-12)
    assert(abs(C_p_min - -0.116094421600922) < 1e-12)
    assert(abs(Cx - 0.142400889444425) < 1e-12)
    assert(abs(Cy) < 1.e-12)
    assert(abs(Cz) < 1e-12)


def test_14_supersonic_half_wing_lower_morino_allow_wake_asym_flow():

    with open("test/input_files/supersonic_half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    input_dict["flow"]["freestream_velocity"] = [100.0, 5.0, 5.0]
    input_dict["geometry"]["wake_model"]["append_wake"] = True

    with open("test/input_files/altered_supersonic_input.json", 'w') as input_handle:
        input_dict = json.dump(input_dict, input_handle, indent=4)

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/altered_supersonic_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.194725933694785) < 1e-12)
    assert(abs(C_p_min - -0.298780589679282) < 1e-12)
    assert(abs(Cx - 0.142781624853592) < 1e-12)
    assert(abs(Cy - 0.000852303050093563) < 1.e-12)
    assert(abs(Cz - 0.892593299837566) < 1e-9)


def test_15_supersonic_half_wing_higher_source_free_zero_aoa_zero_beta():

    with open("test/input_files/supersonic_half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    input_dict["solver"]["formulation"] = 'source-free'
    input_dict["geometry"]["singularity_order"] = "higher"

    with open("test/input_files/altered_supersonic_input.json", 'w') as input_handle:
        input_dict = json.dump(input_dict, input_handle, indent=4)

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/altered_supersonic_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.119055399649912) < 1e-12)
    assert(abs(C_p_min - -0.112090098429673) < 1e-12)
    assert(abs(Cx - 0.142437748460037) < 1e-12)
    assert(abs(Cy) < 1.0e-12)
    assert(abs(Cz) < 1.0e-10)


def test_16_supersonic_half_wing_higher_morino_allow_wake_asym_flow():

    with open("test/input_files/supersonic_half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    input_dict["flow"]["freestream_velocity"] = [100.0, 5.0, 5.0]
    input_dict["geometry"]["wake_model"]["append_wake"] = True
    input_dict["geometry"]["singularity_order"] = "higher"

    with open("test/input_files/altered_supersonic_input.json", 'w') as input_handle:
        input_dict = json.dump(input_dict, input_handle, indent=4)

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/altered_supersonic_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.195559054627881) < 1e-12)
    assert(abs(C_p_min - -0.291045728014071) < 1e-12)
    assert(abs(Cx - 0.142750535947053) < 1e-12)
    assert(abs(Cy - 0.000816008709771361) < 1.e-12)
    assert(abs(Cz - 0.88869372122009) < 1e-9)


def test_17_supersonic_half_wing_source_free_allow_wake_sym_flow():

    with open("test/input_files/supersonic_half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    input_dict["flow"]["freestream_velocity"] = [100.0, 0.0, 5.0]
    input_dict["geometry"]["wake_model"]["append_wake"] = True
    input_dict["solver"]["formulation"] = "source-free"

    with open("test/input_files/altered_supersonic_input.json", 'w') as input_handle:
        input_dict = json.dump(input_dict, input_handle, indent=4)

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/altered_supersonic_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.194950373758414) < 1e-12)
    assert(abs(C_p_min - -0.291642542743292) < 1e-12)
    assert(abs(Cx - 0.142905077607286) < 1e-12)
    assert(abs(Cy) < 1.e-12)
    assert(abs(Cz - 0.893492282700425) < 1e-10)


def test_18_fuselage_subsonic_compressible_iterative():

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/fuselage_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.955084903205978) < 1e-12)
    assert(abs(C_p_min - -0.574040916470699) < 1e-12)
    assert(abs(Cx) < 2.1e-3)
    assert(abs(Cy) < 2.1e-3)
    assert(abs(Cz) < 2.1e-3)


def test_19_supersonic_full_wing_morino_qrup():

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/supersonic_full_wing_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.194950351346633) < 1e-12)
    assert(abs(C_p_min - -0.324945660429385) < 1e-12)
    assert(abs(Cx - 0.0718540012154408) < 1e-12)
    assert(abs(Cy) < 1e-11)
    assert(abs(Cz - 0.429236847680447) < 1e-12)


def test_20_half_wing_lower_morino_asym_inc_flow_off_body():
    # Tests the half wing case with the morino formulation returns consistent off-body results

    # Load original input
    with open("test/input_files/half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Remove old output file
    output_file = "test/results/half_wing_inc_offbody_points.csv"
    try:
        os.remove(output_file)
    except FileNotFoundError:
        pass

    # Alter input
    input_dict["solver"]["formulation"] = "morino"
    input_dict["output"] = {
                                "offbody_points" : {
                                    "points_file" : "test/input_files/root_xz_sample_points.csv",
                                    "output_file" : output_file
                                },
                                "report_file" : "test/results/report.json"
                            }

    # Write altered input
    altered_input_file = "test/input_files/altered_half_wing_input.json"
    with open(altered_input_file, 'w') as altered_input_handle:
        json.dump(input_dict, altered_input_handle, indent=4)

    # Run MachLine
    _ = run_machline("test/input_files/altered_half_wing_input.json", remove_input=True, remove_results=False)

    # Load results
    correct_file = "test/input_files/half_wing_inc_offbody_points_correct.csv"
    results = np.genfromtxt(output_file, skip_header=1, delimiter=',')
    standard = np.genfromtxt(correct_file, skip_header=1, delimiter=',')

    # Update
    if UPDATE_OFFBODY_RESULTS:
        shutil.copyfile(output_file, correct_file)

    # Check
    assert(np.all(np.abs(results - standard) < 1.0e-12))

    # Remove results
    shutil.rmtree("test/results")


def test_21_half_wing_lower_morino_asym_supersonic_flow_off_body():
    # Tests the half wing case with the morino formulation returns consistent off-body results in supersonic flow

    # Load original input
    with open("test/input_files/supersonic_half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Remove old output file
    output_file = "test/results/half_wing_supersonic_offbody_points.csv"
    try:
        os.remove(output_file)
    except FileNotFoundError:
        pass

    # Alter input
    input_dict["flow"]["freestream_velocity"] = [100.0, 5.0, 5.0]
    input_dict["solver"]["formulation"] = "morino"
    input_dict["output"] = {
                                "offbody_points" : {
                                    "points_file" : "test/input_files/root_xz_sample_points.csv",
                                    "output_file" : output_file
                                },
                                "report_file" : "test/results/report.json"
                            }

    # Write altered input
    altered_input_file = "test/input_files/altered_supersonic_half_wing_input.json"
    with open(altered_input_file, 'w') as altered_input_handle:
        json.dump(input_dict, altered_input_handle, indent=4)

    # Run MachLine
    _ = run_machline(altered_input_file, remove_input=True, remove_results=False)

    # Load results
    correct_file = "test/input_files/half_wing_supersonic_offbody_points_correct.csv"
    results = np.genfromtxt(output_file, skip_header=1, delimiter=',')
    standard = np.genfromtxt(correct_file, skip_header=1, delimiter=',')

    # Update
    if UPDATE_OFFBODY_RESULTS:
        shutil.copyfile(output_file, correct_file)

    # Check
    assert(np.all(np.abs(results-standard) < 1.0e-12))

    # Remove results
    shutil.rmtree("test/results")