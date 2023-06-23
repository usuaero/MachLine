import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import write_input_file, run_quad, cases
from studies.paraview_functions import extract_all_data, get_data_column_from_array


RERUN_MACHLINE = False
study_dir = "studies/sphere/"
plot_dir = study_dir + "plots/compressible_pressures/"


def extract_pressures(result_file):

    # Get data
    headers, data = extract_all_data(result_file, which_data='cell')

    # Get locations and pressures
    C_P = get_data_column_from_array(headers, data, 'C_p_ise')
    locs = np.zeros((len(C_P),3))
    locs[:,0] = get_data_column_from_array(headers, data, 'centroid:0') - 1.0
    locs[:,1] = get_data_column_from_array(headers, data, 'centroid:1')
    locs[:,2] = get_data_column_from_array(headers, data, 'centroid:2')

    return locs, C_P


def extract_velocity_magnitudes(result_file):

    # Get data
    headers, data = extract_all_data(result_file, which_data='cell')

    # Get locations and pressures
    u = get_data_column_from_array(headers, data, 'v:0')
    v = get_data_column_from_array(headers, data, 'v:1')
    w = get_data_column_from_array(headers, data, 'v:2')
    V = np.sqrt(u**2 + v**2 + w**2)
    locs = np.zeros((len(u),3))
    locs[:,0] = get_data_column_from_array(headers, data, 'centroid:0') - 1.0
    locs[:,1] = get_data_column_from_array(headers, data, 'centroid:1')
    locs[:,2] = get_data_column_from_array(headers, data, 'centroid:2')

    return locs, V


def get_tamada_pressures():
    # Returns the pressures from Tamada as a function of theta

    # Read in 
    data_file = study_dir + "surface_velocities.csv"
    data = np.genfromtxt(data_file, delimiter=',', skip_header=2)

    # Parse
    theta_04 = np.radians(data[:,0]).flatten()
    v_04 = data[:,1].flatten()
    theta_05 = np.radians(data[:,2]).flatten()
    v_05 = data[:,3].flatten()

    return theta_04, v_04, theta_05, v_05


if __name__=="__main__":

    # Options
    densities = ['ultra_coarse', 'very_coarse', 'coarse', 'medium']
    phis = np.radians([0.0, 30.0, 45.0, 60.0])
    thetas = np.radians([0.0, 30.0, 45.0, 60.0])
    Ms = [0.4, 0.5]

    # Get analytic results
    theta_04, v_04, theta_05, v_05 = get_tamada_pressures()
    theta_anl = [theta_04, theta_05]
    v_anl = [v_04, v_05]

    # Loop
    t = np.zeros((len(Ms), len(densities), len(phis), len(thetas), 4))
    for k, M in enumerate(Ms):
        for l, density in enumerate(densities):

            for i, phi in enumerate(phis):
                for j, theta in enumerate(thetas):

                    # Initialize input
                    result_file = study_dir + "results/sphere_M_{0}_{1}_{2}_{3}.vtk".format(M, round(np.degrees(phi)), round(np.degrees(theta)), density)
                    report_file = study_dir + "reports/sphere_M_{0}_{1}_{2}_{3}.json".format(M, round(np.degrees(phi)), round(np.degrees(theta)), density)

                    V_inf = [np.cos(phi)*np.cos(theta), np.sin(phi)*np.cos(theta), np.sin(theta)]

                    # Lower-order results
                    input_dict ={
                        "flow" : {
                            "freestream_velocity" : V_inf,
                            "gamma" : 1.4,
                            "freestream_mach_number" : M
                        },
                        "geometry" : {
                            "file" : study_dir + "meshes/sphere_{0}.stl".format(density),
                            "spanwise_axis" : "+y"
                        },
                        "solver" : {
                        },
                        "post_processing" : {
                            "pressure_rules" : {
                                "isentropic" : True
                            }
                        },
                        "output" : {
                            "body_file" : result_file,
                            "report_file" : report_file
                        }
                    }

                    # Run
                    input_filename = study_dir + "input.json"
                    write_input_file(input_dict, input_filename)
                    reports = run_quad(input_filename, run=RERUN_MACHLINE)

                    # Get runtimes
                    for m in range(4):
                        t[k,l,i,j,m] = reports[m]["total_runtime"]

                    ## Loop through cases
                    #plt.figure()
                    #for report, case in zip(reports, cases):

                    #    # Get pressures and velocities
                    #    result_file = report["input"]["output"]["body_file"]
                    #    locs, C_P = extract_pressures(result_file)
                    #    _, V = extract_velocity_magnitudes(result_file)

                    #    # Figure out thetas
                    #    theta_space = np.arccos(np.einsum('ij,j->i', locs, V_inf))
                    #    plt.plot(np.degrees(theta_space[np.where(theta_space <= 0.5*np.pi)]), V[np.where(theta_space <= 0.5*np.pi)], 'k.', markersize=1, label='MachLine')
                    #    plt.plot(np.degrees(theta_anl[k]), v_anl[k], 'k-', label='Tamada')
                    #    plt.xlabel('$\\theta [^\circ]$')
                    #    plt.ylabel('$|\mathbf{u}|$')
                    #    if case=='MH':
                    #        plt.legend()
                    #    plt.savefig(plot_dir + "M_{0}_{1}_{2}_{3}_{4}.pdf".format(M, case, round(np.degrees(phi)), round(np.degrees(theta)), density))
                    #    plt.close()

    # Analyze timing
    t_avg = np.average(np.average(np.average(t, axis=0), axis=1), axis=1)
    for i, density in enumerate(densities):
        for j, case in enumerate(cases):
            print(density, case, t_avg[i,j])