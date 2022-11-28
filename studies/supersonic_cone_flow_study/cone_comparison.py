import enum
import json
import subprocess as sp
import numpy as np
import paraview.simple as pvs
import matplotlib.pyplot as plt


def run_comparison(M, grid, half_angle, run_machline=True):
    # Runs the comparison of a cone with the Taylor-MacColl result

    # Parameters
    gamma = 1.4

    # Storage locations
    case_name = "M_{0}_{1}_deg_{2}".format(M, int(half_angle), grid)
    plot_dir = "studies/supersonic_cone_flow_study/plots/"
    mesh_file = "studies/supersonic_cone_flow_study/meshes/cone_{0}_deg_{1}.vtk".format(int(half_angle), grid)
    body_file = "studies/supersonic_cone_flow_study/results/"+case_name+".vtk"
    report_file = "studies/supersonic_cone_flow_study/reports/"+case_name+".json"
    data_file = 'studies/supersonic_cone_flow_study/data/'+case_name+'.csv'

    if run_machline:

        # Declare MachLine input
        input_dict = {
            "flow": {
                "freestream_velocity": [-100.0, 0.0, 0.0],
                "gamma" : gamma,
                "freestream_mach_number" : M
            },
            "geometry": {
                "file": mesh_file,
                "spanwise_axis" : "+y",
                "mirror_about" : "xy",
                "wake_model": {
                    "append_wake" : False,
                },
                "reference": {
                    "area": 4.0
                }
            },
            "solver": {
                "formulation": "morino",
                "control_point_offset": 1.1e-8,
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
                "body_file" : body_file,
                "report_file" : report_file
            }
        }

        # Dump
        input_file = "studies/supersonic_cone_flow_study/cone_input.json"
        with open(input_file, 'w') as input_handle:
            json.dump(input_dict, input_handle, indent=4)

        # Run
        sp.run(["./machline.exe", input_file])
    
    # Load MachLine report file
    try:
        with open(report_file, 'r') as report_handle:
            report = json.load(report_handle)
    except FileNotFoundError:
        return np.nan, np.nan, np.nan, np.nan

    # Read into ParaView
    data_reader = pvs.LegacyVTKReader(registrationName=body_file.replace("dev/results/", ""), FileNames=body_file)

    # Filter cell data to point data
    filter = pvs.CellDatatoPointData(registrationName='Filter', Input=data_reader)
    data_to_process = ['C_p_ise', 'C_p_2nd', 'C_p_sln', 'C_p_lin', 'mu', 'sigma']
    filter.CellDataArraytoprocess = data_to_process

    # Extract and save data
    plot = pvs.PlotOnIntersectionCurves(registrationName='Plot', Input=filter)
    plot.SliceType = 'Plane'
    plot.SliceType.Normal = [0.0, 1.0, 0.0]
    view = pvs.CreateView('XYChartView')
    display = pvs.Show(plot, view, 'XYChartRepresentation')
    view.Update()
    display.XArrayName = 'Points_X'
    view.Update()
    pvs.SaveData(data_file, proxy=plot, PointDataArrays=data_to_process, FieldAssociation='Point Data', Precision=12)

    # Read in data
    data = np.genfromtxt(data_file, delimiter=',', skip_header=1)
    C_p_2nd = data[:,0]
    C_p_ise = data[:,1]
    C_p_lin = data[:,2]
    C_p_sln = data[:,3]
    
    # Plot data from MachLine
    plt.figure()
    plt.plot(data[:,4], C_p_2nd, 'ks', markersize=3, label='Second-Order')
    plt.plot(data[:,4], C_p_ise, 'kv', markersize=3, label='Isentropic')
    plt.plot(data[:,4], C_p_sln, 'ko', markersize=3, label='Slender-Body')
    plt.plot(data[:,4], C_p_lin, 'k^', markersize=3, label='Linear')

    # Format
    plt.xlabel('$x$')
    plt.ylabel('$C_p$')
    plt.legend(title="Pressure Rule", fontsize=6, title_fontsize=6)
    plt.savefig(plot_dir+case_name+".pdf")
    plt.close()

    return C_p_2nd, C_p_ise, C_p_sln, C_p_lin


def get_analytic_data(filename):
    """Reads in the analytic data from the given file."""

    # Get raw data
    raw_data = np.genfromtxt(filename, dtype=str, delimiter=',')

    # Parse data
    Ms = raw_data[2:,1].astype(float)
    thetas =  raw_data[1,2:].astype(float)
    Cps = raw_data[2:,2:].astype(float)

    return Ms, thetas, Cps


if __name__=="__main__":

    # Study parameters
    grids = ["coarse", "medium", "fine"]
    #grids = ["fine"]
    Ms = [1.4, 1.5, 1.7, 2.0, 2.4, 2.8, 3.3, 4.0]
    half_angles = [2.5, 5, 10, 15]

    # Initialize storage
    C_p_2nd_avg = np.zeros((len(grids), len(Ms), len(half_angles)))
    C_p_ise_avg = np.zeros((len(grids), len(Ms), len(half_angles)))
    C_p_sln_avg = np.zeros((len(grids), len(Ms), len(half_angles)))
    C_p_lin_avg = np.zeros((len(grids), len(Ms), len(half_angles)))
    C_p_2nd_s_dev = np.zeros((len(grids), len(Ms), len(half_angles)))
    C_p_ise_s_dev = np.zeros((len(grids), len(Ms), len(half_angles)))
    C_p_sln_s_dev = np.zeros((len(grids), len(Ms), len(half_angles)))
    C_p_lin_s_dev = np.zeros((len(grids), len(Ms), len(half_angles)))

    # Run cases
    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):
            for k, half_angle in enumerate(half_angles):

                if M == 4.0 and half_angle == 15.0:
                    C_p_2nd_avg[i,j,k] = np.nan
                    C_p_ise_avg[i,j,k] = np.nan
                    C_p_sln_avg[i,j,k] = np.nan
                    C_p_lin_avg[i,j,k] = np.nan
                    C_p_2nd_s_dev[i,j,k] = np.nan
                    C_p_ise_s_dev[i,j,k] = np.nan
                    C_p_sln_s_dev[i,j,k] = np.nan
                    C_p_lin_s_dev[i,j,k] = np.nan
                else:
                    C_p_2nd, C_p_ise, C_p_sln, C_p_lin = run_comparison(M, grid, half_angle, run_machline=False)
                    C_p_2nd_avg[i,j,k] = np.average(C_p_2nd).item()
                    C_p_ise_avg[i,j,k] = np.average(C_p_ise).item()
                    C_p_sln_avg[i,j,k] = np.average(C_p_sln).item()
                    C_p_lin_avg[i,j,k] = np.average(C_p_lin).item()
                    C_p_2nd_s_dev[i,j,k] = np.std(C_p_2nd).item()
                    C_p_ise_s_dev[i,j,k] = np.std(C_p_ise).item()
                    C_p_sln_s_dev[i,j,k] = np.std(C_p_sln).item()
                    C_p_lin_s_dev[i,j,k] = np.std(C_p_lin).item()

    # Get analytic data
    Ms_anl, thetas_anl, Cps_anl = get_analytic_data("studies/supersonic_cone_flow_study/Cone Data Zero AoA.csv")

    # Plot overall data
    for k, half_angle in enumerate(half_angles):

        plt.figure()

        # Plots using the most-refined grid
        plt.plot(Ms, C_p_2nd_avg[-1,:,k], mfc='k', marker='s', ls='', mec='k', markersize=3, label='Second-Order')
        plt.plot(Ms, C_p_ise_avg[-1,:,k], mfc='k', marker='v', ls='', mec='k', markersize=3, label='Isentropic')
        plt.plot(Ms, C_p_sln_avg[-1,:,k], mfc='k', marker='o', ls='', mec='k', markersize=3, label='Slender-Body')
        plt.plot(Ms, C_p_lin_avg[-1,:,k], mfc='k', marker='^', ls='', mec='k', markersize=3, label='Linear')

        # Get analytic data for this half angle
        i_theta = np.where(thetas_anl == half_angle)
        plt.plot(Ms_anl[:7], Cps_anl[:7,i_theta].flatten(), 'k', label='Taylor-MacColl', linewidth=1)

        plt.xlabel("$M_\infty$")
        plt.ylabel("$C_p$")
        plt.ylim(bottom=0.0)
        plt.legend(fontsize=6, title_fontsize=6)
        plt.savefig("studies/supersonic_cone_flow_study/plots/C_p_over_M_{0}_deg.pdf".format(half_angle))

    # Determine max standard deviation
    std_max_ise = np.nanmax(np.nanmax(np.nanmax(C_p_ise_s_dev))).item()
    std_max_2nd = np.nanmax(np.nanmax(np.nanmax(C_p_2nd_s_dev))).item()
    std_max_lin = np.nanmax(np.nanmax(np.nanmax(C_p_lin_s_dev))).item()
    std_max_sln = np.nanmax(np.nanmax(np.nanmax(C_p_sln_s_dev))).item()
    print()
    print("Maximum standard deviations:")
    print("    Isentropic: ", std_max_ise)
    print("    Second-order: ", std_max_2nd)
    print("    Linear: ", std_max_lin)
    print("    Slender-body: ", std_max_sln)

    # Plot convergence
    plt.figure()
    for k, half_angle in enumerate(half_angles):
        if k==1:
            plt.plot(Ms, C_p_ise_avg[0,:,k], mfc='w', marker='o', ls='', mec='k', markersize=10, label='Coarse')
            plt.plot(Ms, C_p_ise_avg[1,:,k], mfc='w', marker='o', ls='', mec='k', markersize=7, label='Medium')
            plt.plot(Ms, C_p_ise_avg[2,:,k], mfc='w', marker='o', ls='', mec='k', markersize=4, label='Fine')
        else:
            plt.plot(Ms, C_p_ise_avg[0,:,k], mfc='w', marker='o', ls='', mec='k', markersize=10)
            plt.plot(Ms, C_p_ise_avg[1,:,k], mfc='w', marker='o', ls='', mec='k', markersize=7)
            plt.plot(Ms, C_p_ise_avg[2,:,k], mfc='w', marker='o', ls='', mec='k', markersize=4)
    plt.xlabel("$M_\infty$")
    plt.ylabel("$C_p$")
    plt.ylim(bottom=0.0)
    plt.legend(fontsize=6, title_fontsize=6)
    plt.savefig("studies/supersonic_cone_flow_study/plots/C_p_ise_convergence.pdf")

    plt.figure()
    for k, half_angle in enumerate(half_angles):
        if k==1:
            plt.plot(Ms, C_p_2nd_avg[0,:,k], mfc='w', marker='o', ls='', mec='k', markersize=10, label='Coarse')
            plt.plot(Ms, C_p_2nd_avg[1,:,k], mfc='w', marker='o', ls='', mec='k', markersize=5, label='Medium')
            plt.plot(Ms, C_p_2nd_avg[2,:,k], mfc='w', marker='o', ls='', mec='k', markersize=3, label='Fine')
        else:
            plt.plot(Ms, C_p_2nd_avg[0,:,k], mfc='w', marker='o', ls='', mec='k', markersize=10)
            plt.plot(Ms, C_p_2nd_avg[1,:,k], mfc='w', marker='o', ls='', mec='k', markersize=5)
            plt.plot(Ms, C_p_2nd_avg[2,:,k], mfc='w', marker='o', ls='', mec='k', markersize=3)
    plt.xlabel("$M_\infty$")
    plt.ylabel("$C_p$")
    plt.ylim(bottom=0.0)
    plt.legend(fontsize=6, title_fontsize=6)
    plt.savefig("studies/supersonic_cone_flow_study/plots/C_p_2nd_convergence.pdf")

    plt.figure()
    for k, half_angle in enumerate(half_angles):
        if k==1:
            plt.plot(Ms, C_p_sln_avg[0,:,k], mfc='w', marker='o', ls='', mec='k', markersize=10, label='Coarse')
            plt.plot(Ms, C_p_sln_avg[1,:,k], mfc='w', marker='o', ls='', mec='k', markersize=5, label='Medium')
            plt.plot(Ms, C_p_sln_avg[2,:,k], mfc='w', marker='o', ls='', mec='k', markersize=3, label='Fine')
        else:
            plt.plot(Ms, C_p_sln_avg[0,:,k], mfc='w', marker='o', ls='', mec='k', markersize=10)
            plt.plot(Ms, C_p_sln_avg[1,:,k], mfc='w', marker='o', ls='', mec='k', markersize=5)
            plt.plot(Ms, C_p_sln_avg[2,:,k], mfc='w', marker='o', ls='', mec='k', markersize=3)
    plt.xlabel("$M_\infty$")
    plt.ylabel("$C_p$")
    plt.ylim(bottom=0.0)
    plt.legend(fontsize=6, title_fontsize=6)
    plt.savefig("studies/supersonic_cone_flow_study/plots/C_p_sln_convergence.pdf")

    plt.figure()
    for k, half_angle in enumerate(half_angles):
        if k==1:
            plt.plot(Ms, C_p_lin_avg[0,:,k], mfc='w', marker='o', ls='', mec='k', markersize=10, label='Coarse')
            plt.plot(Ms, C_p_lin_avg[1,:,k], mfc='w', marker='o', ls='', mec='k', markersize=5, label='Medium')
            plt.plot(Ms, C_p_lin_avg[2,:,k], mfc='w', marker='o', ls='', mec='k', markersize=3, label='Fine')
        else:
            plt.plot(Ms, C_p_lin_avg[0,:,k], mfc='w', marker='o', ls='', mec='k', markersize=10)
            plt.plot(Ms, C_p_lin_avg[1,:,k], mfc='w', marker='o', ls='', mec='k', markersize=5)
            plt.plot(Ms, C_p_lin_avg[2,:,k], mfc='w', marker='o', ls='', mec='k', markersize=3)
    plt.xlabel("$M_\infty$")
    plt.ylabel("$C_p$")
    plt.ylim(bottom=0.0)
    plt.legend(fontsize=6, title_fontsize=6)
    plt.savefig("studies/supersonic_cone_flow_study/plots/C_p_lin_convergence.pdf")