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


def test_half_wing_source_free_asym_inc_flow():
    # Tests the half wing case with the source-free formulation returns the consistent result

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/half_wing_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.749203172588596) < 1e-7)
    assert(abs(C_p_min - -1.27770150667225) < 1e-7)
    assert(abs(Cx - -0.393161802358241) < 1e-7)
    assert(abs(Cy - -0.046892410949686) < 1e-7)
    assert(abs(Cz - 20.6218617377561) < 1.1e-7)


def test_half_wing_morino_asym_inc_flow():
    # Tests the half wing case with the source-free formulation returns the consistent result

    # Load original input
    with open("test/half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Alter input
    input_dict["solver"]["formulation"] = "morino"

    # Write altered input
    altered_input_file = "test/altered_half_wing_input.json"
    with open(altered_input_file, 'w') as altered_input_handle:
        json.dump(input_dict, altered_input_handle, indent=4)

    # Run MachLine
    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/altered_half_wing_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.749203884486916) < 1e-7)
    assert(abs(C_p_min - -1.27775688504459) < 1e-7)
    assert(abs(Cx - -0.393164042327084) < 1e-7)
    assert(abs(Cy - -0.0468954335679093) < 1e-7)
    assert(abs(Cz - 20.6220753512639) < 1e-7)


def test_half_wing_morino_zero_aoa_zero_beta_inc():
    # Tests the half wing case at zero angle of attack with the Morino formulation

    # Load original input
    with open("test/half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Alter input
    input_dict["flow"]["freestream_velocity"] = [100.0, 0.0, 0.0]
    input_dict["solver"]["formulation"] = "morino"

    # Write altered input
    altered_input_file = "test/altered_half_wing_input.json"
    with open(altered_input_file, 'w') as altered_input_handle:
        json.dump(input_dict, altered_input_handle, indent=4)

    # Run MachLine
    C_p_max, C_p_min, Cx, Cy, Cz = run_machline(altered_input_file, remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.221658891565984) < 1e-7)
    assert(abs(C_p_min - -0.427587709978507) < 1e-7)
    assert(abs(Cx - 0.301886573119651) < 1e-7)
    assert(abs(Cy) < 1e-12)
    assert(abs(Cz - -1.6953028193344e-07) < 1e-7)


def test_sphere_morino_inc():
    # Tests the sphere case with the Morino formulation returns the consistent result

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/sphere_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.9911419902081444) < 1e-7)
    assert(abs(C_p_min - -1.2378311776012643) < 1e-7)
    assert(abs(Cx - -1.2353806873484363e-07) < 1e-7)
    assert(abs(Cy - 1.557437211145707e-05) < 1e-7)
    assert(abs(Cz - 2.0418661341725652e-06) < 1e-7)


def test_full_half_wing_compare_morino_zero_beta_inc():
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
    print(abs(Cx_full-Cx_half))
    print(abs(Cy_full-Cy_half))
    print(abs(Cz_full-Cz_half))

    assert(abs(C_p_max_full-C_p_max_half)<1e-4)
    assert(abs(C_p_min_full-C_p_min_half)<2e-4)
    assert(abs(Cx_full-Cx_half)<1e-4)
    assert(abs(Cy_full-Cy_half)<1e-4)
    assert(abs(Cz_full-Cz_half)<1e-4)


def test_supersonic_half_wing_morino_zero_aoa_zero_beta():

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/supersonic_half_wing_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.12170587871737629) < 1e-7)
    assert(abs(C_p_min - -0.11606239507721998) < 1e-7)
    assert(abs(Cx - 0.14238137125495076) < 1e-7)
    assert(abs(Cy) < 1e-12)
    assert(abs(Cz) < 1e-7)


def test_supersonic_half_wing_source_free_zero_aoa_zero_beta():

    with open("test/supersonic_half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    input_dict["solver"]["formulation"] = 'source-free'

    with open("test/altered_supersonic_input.json", 'w') as input_handle:
        input_dict = json.dump(input_dict, input_handle, indent=4)

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/altered_supersonic_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.121707118043971) < 1e-7)
    assert(abs(C_p_min - -0.116058405091721) < 1e-7)
    assert(abs(Cx - 0.142379818437667) < 1e-7)
    assert(abs(Cy) < 1.e-12)
    assert(abs(Cz) < 2.1e-7)


def test_supersonic_half_wing_morino_allow_wake_asym_flow():

    with open("test/supersonic_half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    input_dict["flow"]["freestream_velocity"] = [100.0, 5.0, 5.0]
    input_dict["geometry"]["wake_model"]["append_wake"] = True

    with open("test/altered_supersonic_input.json", 'w') as input_handle:
        input_dict = json.dump(input_dict, input_handle, indent=4)

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/altered_supersonic_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.803079382999673) < 1e-7)
    assert(abs(C_p_min - -0.355152163248553) < 1e-7)
    assert(abs(Cx - 0.141371468627056) < 1e-7)
    assert(abs(Cy - 0.000914046488943176) < 1.e-12)
    assert(abs(Cz - 0.910059453750861) < 2.1e-7)


def test_supersonic_half_wing_source_free_allow_wake_sym_flow():

    with open("test/supersonic_half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    input_dict["flow"]["freestream_velocity"] = [100.0, 0.0, 5.0]
    input_dict["geometry"]["wake_model"]["append_wake"] = True

    with open("test/altered_supersonic_input.json", 'w') as input_handle:
        input_dict = json.dump(input_dict, input_handle, indent=4)

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/altered_supersonic_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.194950373758269) < 1e-7)
    assert(abs(C_p_min - -0.301607918132292) < 1e-7)
    assert(abs(Cx - 0.142869043085108) < 1e-7)
    assert(abs(Cy) < 1.e-12)
    assert(abs(Cz - 0.892760941615457) < 2.1e-7)


def test_fuselage_subsonic_compressible_iterative():

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/fuselage_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.955104302502646) < 1e-7)
    assert(abs(C_p_min - -0.56745347779637) < 1e-7)
    assert(abs(Cx) < 2.1e-3)
    assert(abs(Cy) < 2.1e-3)
    assert(abs(Cz) < 2.1e-3)