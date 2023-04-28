import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, cases, line_styles, write_input_file
from studies.paraview_functions import extract_plane_slice, get_data_column_from_array


RERUN_MACHLINE = True
study_dir = "studies/nacelles/"
plot_dir = study_dir + "plots/supersonic/"


def run_study_at_mach(M):
    # Runs the study for the given Mach number

    # options
    grids = ["coarse", "medium", "fine", "ultra_fine"]
    line_colors = ['#BBBBBB', '#999999', '#777777', '#000000']

    # Run MachLine for all results
    report_dicts = []
    for i, grid in enumerate(grids):

        # Files
        mesh_file = study_dir + "meshes/supersonic_{0}.stl".format(grid)
        result_file = study_dir + "results/supersonic_M_{0}_{1}.vtk".format(M, grid)
        control_point_file = study_dir + "results/supersonic_M_{0}_{1}_control_points.vtk".format(M, grid)
        input_file = study_dir + "input.json"
        report_file = study_dir + "reports/supersonic_M_{0}_{1}.json".format(M, grid)

        # Declare input
        input_dict = {
            "flow": {
                "freestream_velocity": [1.0, 0.0, 0.0],
                "gamma" : 1.4,
                "freestream_mach_number" : M
            },
            "geometry": {
                "file": mesh_file,
                "spanwise_axis" : "+y",
                "max_continuity_angle" : 5.0,
                "wake_model": {
                    "append_wake" : False,
                },
                "reference": {
                }
            },
            "solver": {
            },
            "post_processing" : {
                "pressure_rules" : {
                    "isentropic" : True
                }
            },
            "output" : {
                "body_file" : result_file,
                "control_point_file" : control_point_file,
                "report_file" : report_file
            }
        }

        # Write input
        write_input_file(input_dict, input_file)

        # Run
        reports = run_quad(input_file, run=RERUN_MACHLINE)
        report_dicts.append(list(reports))

    # Plot

    # Loop through cases
    for j, case in enumerate(cases):

        # Initialize figures
        inner_fig = plt.figure(1)
        inner_ax = inner_fig.subplots()
        outer_fig = plt.figure(2)
        outer_ax = outer_fig.subplots()

        # Loop through grid densities
        for i, grid in enumerate(grids):

            report = report_dicts[i][j]

            # Get result file
            result_file = report["input"]["output"]["body_file"]

            # Load slice
            headers, data = extract_plane_slice(result_file, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], which_data='cell')

            # Get data
            x = get_data_column_from_array(headers, data, 'centroid:0')
            z = get_data_column_from_array(headers, data, 'centroid:2')
            C_P = get_data_column_from_array(headers, data, 'C_p_ise')

            # Splice in sides
            N = len(x)
            X = np.zeros(N)
            X[::2] = x[:N//2]
            X[1::2] = x[N//2:]
            x = X
            CP = np.zeros(N)
            CP[::2] = C_P[:N//2]
            CP[1::2] = C_P[N//2:]
            C_P = CP

            # Get location of front
            front_ind = np.argmin(x) + 2

            # Plot inside pressures
            inner_ax.plot(x[:front_ind], C_P[:front_ind], color=line_colors[i])

            # Plot outside pressures
            outer_ax.plot(x[front_ind:], C_P[front_ind:], color=line_colors[i])
            #plt.figure(1)
            #plt.plot(x)
            #plt.savefig('test.pdf')
            #plt.close(1)

        # Format plots
        inner_ax.set_xlabel("$x$")
        inner_ax.set_ylabel("$C_{P_{ise}}$")
        inner_fig.savefig(plot_dir + "inside_pressures_M_{0}_{1}.pdf".format(M, case))
        plt.close(inner_fig)

        outer_ax.set_xlabel("$x$")
        outer_ax.set_ylabel("$C_{P_{ise}}$")
        outer_fig.savefig(plot_dir + "outside_pressures_M_{0}_{1}.pdf".format(M, case))
        plt.close(outer_fig)


if __name__=="__main__":

    # Mach numbers
    Ms = [1.5, 1.9, 2.0, 2.1, 2.5]

    for M in Ms:
        run_study_at_mach(M)