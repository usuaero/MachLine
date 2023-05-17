import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, write_input_file, cases, line_styles

RERUN_MACHLINE = True
study_dir = "studies/supersonic_spindle/"
mesh_dir = study_dir + "meshes/"
report_dir = study_dir + "reports/"

if __name__=="__main__":

    # Double-wedge
    input_dict = {
        "flow": {
            "freestream_velocity": [ 694.4475790151479, 0.0, 0.0 ],
            "gamma": 1.4,
            "freestream_mach_number": 2.0
        },
        "geometry": {
            "file": mesh_dir + "ehlers_spindle_fine.vtk",
            "spanwise_axis": "+y",
            "wake_model": {
                "append_wake": False
            },
            "reference": {
                "area": 4.0
            }
        },
        "solver": {
            "formulation": "source-free",
            "control_point_offset": 1.1e-08,
            "control_point_offset_type" : 'direct'
        },
        "post_processing": {
            "pressure_rules": {
                "isentropic": True
            }
        },
        "output": {
            "report_file" : report_dir + "report.json"
        }
    }

    # Set up offsets
    offsets = np.logspace(-11.999, 0.001, 13)
    Cx = np.zeros((4,13))
    Cy = np.zeros((4,13))
    Cz = np.zeros((4,13))

    for i, offset in enumerate(offsets):

        # Set offset
        input_dict["solver"]["control_point_offset"] = offset
        input_file = study_dir + "input.json"
        write_input_file(input_dict, input_file)

        # Run
        reports = run_quad(input_file, run=RERUN_MACHLINE)
        
        # Load results
        for j, report in enumerate(reports):
            if report is None:
                Cx[j,i] = np.nan
                Cy[j,i] = np.nan
                Cz[j,i] = np.nan
            else:
                Cx[j,i] = report["total_forces"]["Cx"]
                Cy[j,i] = report["total_forces"]["Cy"]
                Cz[j,i] = report["total_forces"]["Cz"]

    # Plot y and z coefficients
    plt.figure()
    for j, (case, line_style) in enumerate(zip(cases, line_styles)):
        plt.plot(offsets, np.abs(Cy[j]), line_style, label=case)
        plt.plot(offsets, np.abs(Cz[j]), line_style)
    plt.xlabel('$k_1$')
    plt.ylabel('Off-Axis Force Coefficients')
    plt.xscale('log')
    plt.yscale('log')
    #plt.legend()
    plt.show()

    # Plot x coefficient
    plt.figure()
    for j, (case, line_style) in enumerate(zip(cases, line_styles)):
        plt.plot(offsets, np.abs(Cx[j]), line_style, label=case)
    plt.xlabel('$k_1$')
    plt.ylabel('Axial Force Coefficient')
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend()
    plt.show()