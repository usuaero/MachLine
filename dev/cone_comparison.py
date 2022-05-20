import json
import subprocess as sp
import flow54 as fl
import numpy as np
import paraview.simple as pvs
import matplotlib.pyplot as plt


def run_comparison(M, grid, half_angle, run_machline=True):
    # Runs the comparison of a cone with the Taylor-MacColl result

    # Parameters
    gamma = 1.4

    # Storage locations
    case_name = "M_{0}_{1}_deg_{2}".format(M, int(half_angle), grid)
    plot_dir = "dev/results/cone_comparison/plots/"
    body_file = "dev/results/cone_comparison/results/"+case_name+".vtk"
    report_file = "dev/results/cone_comparison/reports/"+case_name+".json"
    data_file = 'dev/results/cone_comparison/data/'+case_name+'.csv'

    if run_machline:

        # Declare MachLine input
        input_dict = {
            "flow": {
                "freestream_velocity": [100.0, 0.0, 0.0],
                "gamma" : gamma,
                "freestream_mach_number" : M
            },
            "geometry": {
                "file": "dev/results/cone_comparison/meshes/cone_{0}_deg_{1}.vtk".format(int(half_angle), grid),
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
                "control_point_file" : body_file.replace('.vtk', '_control_points.vtk'),
                "report_file" : report_file
            }
        }

        # Dump
        input_file = "dev/results/cone_comparison/cone_input.json"
        with open(input_file, 'w') as input_handle:
            json.dump(input_dict, input_handle, indent=4)

        # Run
        sp.run(["./machline.exe", input_file])
    
    # Load MachLine report file
    with open(report_file, 'r') as report_handle:
        report = json.load(report_handle)

    # Read into ParaView
    data_reader = pvs.LegacyVTKReader(registrationName=body_file.replace("dev/results/", ""), FileNames=body_file)

    # Filter cell data to point data
    filter = pvs.CellDatatoPointData(registrationName='Filter', Input=data_reader)
    data_to_process = ['C_p_ise', 'C_p_2nd', 'C_p_sln', 'C_p_lin', 'mu', 'sigma']
    filter.CellDataArraytoprocess = data_to_process

    # Extract and save data
    plot = pvs.PlotOnIntersectionCurves(registrationName='Plot', Input=filter)
    plot.SliceType = 'Plane'
    plot.SliceType.Normal = [0.0, 0.0, 1.0]
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
    C_p_sln = data[:,2]
    C_p_lin = data[:,3]
    
    # Plot data from MachLine
    plt.figure()
    plt.plot(data[:,4], C_p_2nd, 'ks', markersize=3, label='2nd')
    plt.plot(data[:,4], C_p_ise, 'kv', markersize=3, label='Ise.')
    plt.plot(data[:,4], C_p_sln, 'ko', markersize=3, label='Slnd.')
    plt.plot(data[:,4], C_p_lin, 'k^', markersize=3, label='Lin.')

    # Format
    plt.xlabel('$x$')
    plt.ylabel('$C_p$')
    plt.legend()
    plt.savefig(plot_dir+case_name+".pdf")
    plt.close()

    return C_p_2nd, C_p_ise, C_p_sln, C_p_lin


if __name__=="__main__":

    grids = ["coarse", "medium", "fine"]
    N_verts = [402, 1090, 3770, 13266]
    Ms = [1.5, 2.0, 3.0, 5.0]
    half_angles = [1, 5, 10, 15]

    C_p_2nd = np.zeros((len(grids), len(Ms), len(half_angles)))
    C_p_ise = np.zeros((len(grids), len(Ms), len(half_angles)))
    C_p_sln = np.zeros((len(grids), len(Ms), len(half_angles)))
    C_p_lin = np.zeros((len(grids), len(Ms), len(half_angles)))

    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):
            for k, half_angle in enumerate(half_angles):

                C_p_2nd[i,j,k], C_p_ise[i,j,k], C_p_sln[i,j,k], C_p_lin[i,j,k] = run_comparison(M, grid, half_angle)