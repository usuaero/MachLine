import json
import subprocess as sp
import flow54 as fl
import numpy as np
import paraview.simple as pvs
import matplotlib.pyplot as plt


if __name__=="__main__":

    # Parameters
    M = 1.5
    R_G = 287.058
    gamma = 1.4
    T_inf = 300.0
    c_inf = np.sqrt(gamma*R_G*T_inf)
    alpha = 0.0
    grid = 9800 # 9800

    # Declare MachLine input
    body_file = "dev/results/sears_haack_{0}.vtk".format(grid)
    input_dict = {
        "flow": {
            "freestream_velocity": [M*c_inf*np.cos(np.radians(alpha)), 0.0, M*c_inf*np.sin(np.radians(alpha))],
            "gamma" : gamma,
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": "dev/meshes/sears_haack_{0}.tri".format(grid),
            "spanwise_axis" : "+y",
            "wake_model": {
                "wake_present" : False,
            },
            "reference": {
                "area": 1.0
            }
        },
        "solver": {
            "formulation": "morino",
            "control_point_offset": 1.1e-05
        },
        "post_processing" : {
            "pressure_rules" : {
                "second-order" : True,
                "isentropic" : True
            }
        },
        "output" : {
            "body_file" : body_file,
            "control_point_file" : "dev/results/sears_haack_{0}_control_points.vtk".format(grid)
        }
    }

    # Dump
    input_file = "dev/sears_haack_input.json"
    with open(input_file, 'w') as input_handle:
        json.dump(input_dict, input_handle, indent=4)

    # Run
    sp.run(["./machline.exe", input_file])
    
    # Read into ParaView
    data_reader = pvs.LegacyVTKReader(registrationName=body_file.replace("dev/results/", ""), FileNames=body_file)

    # Filter cell data to point data
    filter = pvs.CellDatatoPointData(registrationName='Filter', Input=data_reader)
    data_to_process = ['C_p_ise', 'C_p_2nd', 'mu', 'sigma']
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
    save_loc = 'dev/results/sears_haack_{0}_data.csv'.format(grid)
    pvs.SaveData(save_loc, proxy=plot, PointDataArrays=data_to_process, FieldAssociation='Point Data', Precision=12)

    # Read in data
    data = np.genfromtxt(save_loc, delimiter=',', skip_header=1)
    
    # Plot data from MachLine
    plt.figure()
    plt.plot(data[:,2], data[:,0], 'ks-', label='MachLine 2nd-order')
    plt.plot(data[:,2], data[:,1], 'kv-', label='MachLine Isentropic')

    # Format
    plt.xlabel('$x$')
    plt.ylabel('$C_p$')
    plt.legend()
    plt.show()