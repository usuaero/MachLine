import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, write_input_file, cases, line_styles

RERUN_MACHLINE = False
study_dir = "studies/supersonic_agard_b_wing_body/"
plot_dir = study_dir + "plots/"


def run_quad_for_mach_and_mesh(M, grid):
    """Runs a case quad for the given Mach number and mesh density."""

    # Storage locations
    case_name = "M_{0}_{1}".format(round(M, 3), grid)
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
                "area": 1.0
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

    # Experimental data from combined study
    exp_data = np.genfromtxt(study_dir + "wind_tunnel_data.csv", delimiter=',', skip_header=2)
    M_exp = exp_data[:,0]
    CD_exp = exp_data[:,1]
    upper_unc = exp_data[:,3] - CD_exp
    lower_unc = CD_exp - exp_data[:,5]

    # Extended experimental data
    M_ext = [1.26, 1.4, 1.6]
    CD_ext = [0.027026848265697, 0.026054637036194, 0.02480864124969]
    CP_b_ext = [-0.246781115879828, -0.215076335877863, -0.195737342804789]

    # Reference area
    D = 0.64516 # Model diameter
    S_ref = 4.0*np.sqrt(3.0)*D**2

    # Base pressure calcs
    S_base = 0.25*np.pi*D**2

    # Parameters
    grids = ["coarse", "medium", "fine"]
    #Ms = np.linspace(1.1, 2.0, 10)
    Ms = list(M_exp[1:]) + M_ext
    CP_b = np.zeros_like(Ms)
    CP_b[-3:] = CP_b_ext
    CD = np.zeros((len(grids), len(Ms), 4))

    # Loop
    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):

            _,_,C_F = run_quad_for_mach_and_mesh(M, grid)
            CD[i,j,:] = -C_F[:,2]/S_ref

    # Plot
    for i, (case, line_style) in enumerate(zip(cases, line_styles)):
        plt.figure()
        for j, grid in enumerate(grids):
            plt.plot(Ms, CD[j,:,i] - CP_b*S_base/S_ref, 'ko', mfc='none', markersize=10-3*j, label="ML "+grid.title())
        plt.errorbar(M_exp, CD_exp, yerr=np.array([lower_unc, upper_unc]), fmt='ks', label='Damljanovic et al. (2020)', markersize=3)
        plt.plot(M_ext, CD_ext, 'kv', label='Damljanovic et al. (2006)')
        plt.xlabel('$M_\infty$')
        plt.ylabel('$C_x$')
        #plt.legend()
        plt.savefig(plot_dir+"c_z_over_M_{0}.pdf".format(case))
        plt.savefig(plot_dir+"c_z_over_M_{0}.svg".format(case))
        plt.close()