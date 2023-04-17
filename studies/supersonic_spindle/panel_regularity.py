import numpy as np
import matplotlib.pyplot as plt

from dev.helper_scripts.geometry_creator import generate_random_spindle
from studies.paraview_functions import extract_all_data, get_data_column_from_array
from studies.case_running_functions import run_quad, write_input_file, cases


RERUN_MACHLINE = True
REGENERATE_MESHES = True

# Parameters
study_dir = "studies/supersonic_spindle/"
input_file = study_dir + "spindle_input.json"
data_file = study_dir + "results/random_spindle_data.csv"
analytic_data_file = study_dir + "method_of_char_from_ehlers.csv"
plot_dir = study_dir + "plots/regularity/"
M = np.sqrt(2)


def run_random_spindle(N, cpo=1.1e-6, regenerate_mesh=False):
    """Runs a random spindle with N vertices and a control point offset of cpo through MachLine and plots the pressures."""

    # Name results file for this number of vertices
    mesh_file = study_dir + "meshes/random_spindle_{0}.vtk".format(N)
    result_file = study_dir + "results/random_spindle_{0}_{1}.vtk".format(N, round(np.log10(cpo)))
    report_file = study_dir + "reports/random_spindle_{0}_{1}.json".format(N, round(np.log10(cpo)))

    # Generate mesh
    if regenerate_mesh:
        def r_of_x(x):
            return 0.2*x*(1.0-x)
        generate_random_spindle(mesh_file, N, 1.0, r_of_x)

    # Write input
    input_dict = {
        "flow": {
            "freestream_velocity": [1.0, 0.0, 0.0],
            "gamma" : 1.4,
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "wake_model" : {
                "wake_present" : False
            }
        },
        "solver" : {
            "matrix_solver" : "GMRES",
            "control_point_offset" : cpo,
            "control_point_offset_type" : "direct"
        },
        "post-processing" : {
            "pressure_rules"  : {
                "isentropic" : True,
                "second-order" : True,
                "slender-body" : True,
                "linear" : True
            }
        },
        "output": {
            "body_file" : result_file,
            "control_point_file" : result_file.replace(".vtk", "_control_points.vtk"),
            "report_file" : report_file
        }
    }
    write_input_file(input_dict, input_file)

    # Run
    reports = run_quad(input_file, run=RERUN_MACHLINE, delete_input=True)

    for i, (case, report) in enumerate(zip(cases, reports)):

        if report is None:
            continue

        if report["solver_results"]["solver_status_code"] != 0:
            continue

        # Get result file
        quad_result_file = report["input"]["output"]["body_file"]

        # Extract data from MachLine
        headers, data = extract_all_data(quad_result_file, which_data='cell')
        x_ML = get_data_column_from_array(headers, data, 'centroid:0')
        C_p_ML = get_data_column_from_array(headers, data, 'C_p_ise')

        # Extract analytic data
        ML_data = np.genfromtxt(analytic_data_file, delimiter=',', skip_header=2)
        C_p_anl = ML_data[:,1]
        x_anl = ML_data[:,0]

        # Plot
        plt.figure()
        plt.plot(x_ML, C_p_ML, 'k.', markersize=3, label='MachLine')
        plt.plot(x_anl, C_p_anl, 'k', label='Meth. of Char.')
        plt.xlabel("$x$")
        plt.ylabel("$C_P$")
        plt.legend(fontsize=6, title_fontsize=6)
        plt.savefig(plot_dir + 'M_{0}_{1}_{2}_{3}.pdf'.format(round(M, 2), round(np.log10(cpo)), N, case))
        plt.savefig(plot_dir + 'M_{0}_{1}_{2}_{3}.svg'.format(round(M, 2), round(np.log10(cpo)), N, case))
        plt.close()


if __name__=="__main__":

    # Numbers of vertices
    Ns = [50, 100, 200, 400, 800]
    
    # Control point offsets
    cpos = [1.1e-4, 1.1e-6, 1.1e-8, 1.1e-10]

    for N in Ns:
        for cpo in cpos:
            run_random_spindle(N, cpo, regenerate_mesh=REGENERATE_MESHES)