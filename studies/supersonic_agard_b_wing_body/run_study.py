import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, write_input_file, cases, line_styles

RERUN_MACHLINE = True
study_dir = "studies/supersonic_agard_b_wing_body/"
plot_dir = study_dir + "plots/"


def run_quad_for_mach_and_mesh(M, grid):
    """Runs a case quad for the given Mach number and mesh density."""

    # Storage locations
    case_name = "M_{0}_{1}".format(M, grid)
    mesh_file = study_dir + "meshes/agard_b_{0}.vtk".format(grid)
    results_file = study_dir + "results/"+case_name+".vtk"
    control_point_file = study_dir + "results/"+case_name+"_control_points.vtk"
    report_file = study_dir + "reports/"+case_name+".json"

    # Write out input file

    # Declare MachLine input
    input_dict = {
        "flow": {
            "freestream_velocity": [0.0, 0.0, -1.0],
            "gamma" : 1.4,
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+x",
            "mirror_about" : "yz",
            "wake_model": {
            },
            "reference": {
                "area": 0.0929
            }
        },
        "solver": {
            "formulation": "morino",
            "run_checks" : True
        },
        "post_processing" : {
            "pressure_rules" : {
                "second-order" : True,
                "isentropic" : True,
                "slender-body" : True,
                "linear" : True
            }
        },
        "output" : {
            "body_file" : results_file,
            "control_point_file" : control_point_file,
            "report_file" : report_file
        }
    }

    # Dump
    input_file = study_dir + "input.json"
    write_input_file(input_dict, input_file)

    # Run quad
    reports = run_quad(input_file, run=RERUN_MACHLINE)

    # Pull out forces
    C_F = np.zeros((4,3))
    for i, report in enumerate(reports):
        try:
            C_F[i,0] = report["total_forces"]["Cx"]
            C_F[i,1] = report["total_forces"]["Cy"]
            C_F[i,2] = report["total_forces"]["Cz"]
        except KeyError:
            C_F[i,0] = np.nan
            C_F[i,1] = np.nan
            C_F[i,2] = np.nan

    # Get system dimension and average characteristic length
    N_sys = reports[0]["solver_results"]["system_dimension"]
    l_avg = reports[0]["mesh_info"]["average_characteristic_length"]

    return N_sys, l_avg, C_F



if __name__=="__main__":

    # Experimental data
    axial_force_FS = 1130 # N
    axial_force_unc = 0.002*axial_force_FS

    # Parameters
    grids = ["coarse", "medium", "fine"]
    Ms = [1.01, 1.26, 1.4, 1.6]
    Cz = np.zeros((len(grids), len(Ms), 4))

    # Loop
    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):

            _,_,C_F = run_quad_for_mach_and_mesh(M, grid)
            Cz[i,j,:] = C_F[:,2]

    # Plot
    for i, (case, line_style) in enumerate(zip(cases, line_styles)):
        plt.figure()
        for j, grid in enumerate(grids):
            plt.plot(Ms, -Cz[j,:,i], 'ko', mfc='white', markersize=10-3*j)
        plt.xlabel('$M_\infty$')
        plt.ylabel('$C_x$')
        plt.savefig(plot_dir+"c_z_over_M_{0}.pdf".format(case))
        plt.close()