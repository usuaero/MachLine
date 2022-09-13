import os
import json
import shutil
import subprocess as sp


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


def test_half_wing_source_free_asym_inc_flow():
    # Tests the half wing case with the source-free formulation returns the consistent result

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/half_wing_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.749203172758166) < 1e-7)
    assert(abs(C_p_min - -1.27767126744922) < 1e-7)
    assert(abs(Cx - -0.393161867046978) < 1e-7)
    assert(abs(Cy - -0.0468922101835385) < 1e-7)
    assert(abs(Cz - 20.6218616569872) < 1.1e-7)


def test_half_wing_morino_asym_inc_flow():
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

    assert(abs(C_p_max - 0.749203890042168) < 1e-7)
    assert(abs(C_p_min - -1.27772665490914) < 1e-7)
    assert(abs(Cx - -0.393164107436776) < 1e-7)
    assert(abs(Cy - -0.0468952332211119) < 1e-7)
    assert(abs(Cz - 20.6220753189898) < 1e-7)


def test_half_wing_morino_zero_aoa_zero_beta_inc():
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

    assert(abs(C_p_max - 0.221659165930464) < 1e-7)
    assert(abs(C_p_min - -0.42758575480925) < 1e-7)
    assert(abs(Cx - 0.30188654863457) < 1e-7)
    assert(abs(Cy) < 1e-12)
    assert(abs(Cz) < 2.2e-7)


def test_sphere_morino_inc_quad_doublets_iterative():
    # Tests the sphere case with the Morino formulation using quadratic doublets and an iterative matrix solver

    # Load original input
    with open("test/input_files/sphere_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Alter input
    input_dict["solver"]["matrix_solver"] = "BJAC"
    input_dict["solver"]["relaxation"] = 0.9
    input_dict["geometry"]["singularity_order"] = {"doublet" : 2}

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

    assert(abs(C_p_max - 0.996050603246009) < 1e-7)
    assert(abs(C_p_min - -1.24733544905577) < 1e-7)
    assert(abs(Cx) < 1e-6)
    assert(abs(Cy) < 1e-6)
    assert(abs(Cz) < 1e-6)


def test_sphere_morino_inc():
    # Tests the sphere case with the Morino formulation returns the consistent result

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/sphere_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.991141910232862) < 1e-7)
    assert(abs(C_p_min - -1.2378315549138) < 1e-7)
    assert(abs(Cx) < 1e-4)
    assert(abs(Cy) < 1e-4)
    assert(abs(Cz) < 1e-4)


def test_full_half_wing_compare_morino_zero_beta_inc():
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


def test_subsonic_comp_pressure_corrections():

    # Tests the half wing case with the source-free formulation returns the consistent result

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/compressible_half_wing_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.817136935804372) < 1e-7)
    assert(abs(C_p_min - -1.09340259165792) < 1e-7)
    assert(abs(Cx - -0.507724739979628) < 1e-7)
    assert(abs(Cy - 0.0101999369033425) < 1e-7)
    assert(abs(Cz - 26.2478893799577) < 1.1e-7)


def test_supersonic_half_wing_morino_zero_aoa_zero_beta():

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/supersonic_half_wing_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.121696594166226) < 1e-7)
    assert(abs(C_p_min - -0.116093853140189) < 1e-7)
    assert(abs(Cx - 0.142400846628573) < 1e-7)
    assert(abs(Cy) < 1e-12)
    assert(abs(Cz) < 1e-6)


def test_supersonic_half_wing_source_free_zero_aoa_zero_beta():

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

    assert(abs(C_p_max - 0.121696594161608) < 1e-7)
    assert(abs(C_p_min - -0.116093852996235) < 1e-7)
    assert(abs(Cx - 0.142400846623555) < 1e-7)
    assert(abs(Cy) < 1.e-12)
    assert(abs(Cz) < 1e-6)


def test_supersonic_half_wing_morino_allow_wake_asym_flow():

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

    assert(abs(C_p_max - 0.194725933709647) < 1e-7)
    assert(abs(C_p_min - -0.298886927641986) < 1e-7)
    assert(abs(Cx - 0.142766060995649) < 1e-7)
    assert(abs(Cy - 0.000849118583451407) < 1.e-11)
    assert(abs(Cz - 0.891813872983368) < 2.1e-7)


def test_supersonic_half_wing_source_free_allow_wake_sym_flow():

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

    assert(abs(C_p_max - 0.19495037375827) < 1e-7)
    assert(abs(C_p_min - -0.29150016335519) < 1e-7)
    assert(abs(Cx - 0.142893375049277) < 1e-7)
    assert(abs(Cy) < 1.e-12)
    assert(abs(Cz - 0.892752615765537) < 2.1e-7)


def test_fuselage_subsonic_compressible_iterative():

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/fuselage_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.955104302502646) < 1e-7)
    assert(abs(C_p_min - -0.567453182856627) < 1e-7)
    assert(abs(Cx) < 2.1e-3)
    assert(abs(Cy) < 2.1e-3)
    assert(abs(Cz) < 2.1e-3)


def test_subsonic_quad_doublets_naca_wing_asym():
    # Tests subsonic quadartic doublets on a half wing in asymmetric flow

    # Load original input
    with open("test/input_files/half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    # Alter input
    input_dict["flow"]["freestream_velocity"] = [100.0, 10.0, 10.0]
    input_dict["solver"]["formulation"] = "morino"
    input_dict["geometry"]["singularity_order"] = { "doublet" : 2 }
    input_dict["geometry"]["file"] = "test/meshes/naca_0010_AR_10_half_coarse.stl"

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

    assert(abs(C_p_max - 0.994093723964344) < 1e-7)
    assert(abs(C_p_min - -9.75794552659239) < 1e-6)
    assert(abs(Cx - -0.228460475171211) < 1e-7)
    assert(abs(Cy - 0.00272511910436683) < 1e-7)
    assert(abs(Cz - 3.84392498845133) < 1e-7)


def test_supersonic_half_wing_morino_quad_doublets_lin_sources_asym_flow():

    with open("test/input_files/supersonic_half_wing_input.json", 'r') as input_handle:
        input_dict = json.load(input_handle)

    input_dict["geometry"]["singularity_order"] = {
        "source" : 1,
        "doublet" : 2
    }
    input_dict["flow"]["freestream_velocity"] = [100.0, 5.0, 5.0]
    input_dict["output"]["body_file"] = "test/results/supersonic_half_wing_right.vtk"
    input_dict["output"]["mirrored_body_file"] = "test/results/supersonic_half_wing_left.vtk"

    with open("test/input_files/altered_supersonic_input.json", 'w') as input_handle:
        input_dict = json.dump(input_dict, input_handle, indent=4)

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/altered_supersonic_input.json", remove_input=True)

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.313596680697127) < 1e-7)
    assert(abs(C_p_min - -0.357142857142857) < 1e-7)
    assert(abs(Cx - 0.129560418577158) < 1e-7)
    assert(abs(Cy - 0.00108973650470752) < 1e-7)
    assert(abs(Cz - 0.826402787381577) < 2e-7)


def test_supersonic_full_wing_morino_qrup():

    C_p_max, C_p_min, Cx, Cy, Cz = run_machline("test/input_files/supersonic_full_wing_input.json")

    print(C_p_max)
    print(C_p_min)
    print(Cx)
    print(Cy)
    print(Cz)

    assert(abs(C_p_max - 0.194950351346633) < 1e-7)
    assert(abs(C_p_min - -0.32494014106503) < 1e-7)
    assert(abs(Cx - 0.0718531668938911) < 1e-7)
    assert(abs(Cy) < 1e-12)
    assert(abs(Cz - 0.429190483258204) < 1e-7)