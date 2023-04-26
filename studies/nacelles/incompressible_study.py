import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, cases, line_styles, write_input_file
from studies.paraview_functions import extract_plane_slice, get_data_column_from_array


RERUN_MACHLINE = False
study_dir = "studies/nacelles/"
plot_dir = study_dir + "plots/incompressible/"


def run_study():
    # Runs the study

    # options
    grids = ["coarse", "medium", "fine"]
    line_colors = ['#BBBBBB', '#888888', '#000000']

    # Run MachLine for all results
    report_dicts = []
    for i, grid_circ in enumerate(grids):
        report_dicts.append([])
        for j, grid_ax in enumerate(grids):

            # Files
            mesh_file = study_dir + "meshes/NACA_0010_{0}_{1}.stl".format(grid_circ, grid_ax)
            result_file = study_dir + "results/NACA_0010_{0}_{1}.vtk".format(grid_circ, grid_ax)
            control_point_file = study_dir + "results/NACA_0010_{0}_{1}_control_points.vtk".format(grid_circ, grid_ax)
            input_file = study_dir + "input.json"
            report_file = study_dir + "reports/NACA_0010_{0}_{1}.json".format(grid_circ, grid_ax)

            # Declare input
            input_dict = {
                "flow": {
                    "freestream_velocity": [1.0, 0.0, 0.0]
                },
                "geometry": {
                    "file": mesh_file,
                    "spanwise_axis" : "+y",
                    "wake_model": {
                        "wake_present" : True,
                        "append_wake" : False
                    },
                    "reference": {
                    }
                },
                "solver": {
                },
                "post_processing" : {
                    "pressure_rules" : {
                        "incompressible" : True
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
            report_dicts[i].append(list(reports))

    # Plot

    # Loop through cases
    for k, case in enumerate(cases):

        # Initialize figures
        axial_fig = plt.figure(1)
        axial_ax = axial_fig.subplots()
        circu_fig = plt.figure(2)
        circu_ax = circu_fig.subplots()

        # Loop through grid densities
        for i, grid in enumerate(grids):

            # Get results for varied circumferential grid
            report = report_dicts[i][-1][k]

            # Get result file
            result_file = report["input"]["output"]["body_file"]

            # Load slice
            headers, data = extract_plane_slice(result_file, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], which_data='cell')

            # Get data
            x = get_data_column_from_array(headers, data, 'centroid:0')
            C_P = get_data_column_from_array(headers, data, 'C_p_inc')

            # Plot varied circumferential pressures
            circu_ax.plot(x, C_P, color=line_colors[i])

            # Get results for varied axial grid
            report = report_dicts[-1][i][k]

            # Get result file
            result_file = report["input"]["output"]["body_file"]

            # Load slice
            headers, data = extract_plane_slice(result_file, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], which_data='cell')

            # Get data
            x = get_data_column_from_array(headers, data, 'centroid:0')
            C_P = get_data_column_from_array(headers, data, 'C_p_inc')

            # Plot varied axial pressures
            axial_ax.plot(x, C_P, color=line_colors[i])
            #plt.figure(1)
            #plt.plot(x)
            #plt.savefig('test.pdf')
            #plt.close(1)

        # Format plots
        axial_ax.set_xlabel("$x$")
        axial_ax.set_ylabel("$C_{P_{inc}}$")
        axial_fig.savefig(plot_dir + "axial_pressures_{0}.pdf".format(case))
        plt.close(axial_fig)

        circu_ax.set_xlabel("$x$")
        circu_ax.set_ylabel("$C_{P_{inc}}$")
        circu_fig.savefig(plot_dir + "circumferential_pressures_{0}.pdf".format(case))
        plt.close(circu_fig)


if __name__=="__main__":

    run_study()