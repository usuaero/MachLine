import json
import subprocess as sp
import flow54 as fl
import numpy as np
import paraview.simple as pvs
import matplotlib.pyplot as plt


if __name__=="__main__":

    # Parameters
    M = 2.0
    R_G = 287.058
    gamma = 1.4
    T_inf = 300.0
    c_inf = np.sqrt(gamma*R_G*T_inf)
    alpha = 0.0
    rho = 1.225
    p_inf = 1.0e5
    grid = "ultra_fine" # coarse, medium, fine, ultra_fine

    # Declare MachLine input
    body_file = "dev/results/diamond_5_deg_full_{0}.vtk".format(grid)
    input_dict = {
        "flow": {
            "freestream_velocity": [M*c_inf*np.cos(np.radians(alpha)), 0.0, M*c_inf*np.sin(np.radians(alpha))],
            "gamma" : gamma,
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": "dev/meshes/diamond_5_deg_full_{0}.stl".format(grid),
            "spanwise_axis" : "+y",
            "wake_model": {
                "append_wake" : False,
            },
            "reference": {
                "area": 4.0
            }
        },
        "solver": {
            "formulation": "morino",
            #"formulation": "source-free",
            "control_point_offset": 1.1e-5
        },
        "post_processing" : {
            "pressure_rules" : {
                "second-order" : True,
                "isentropic" : True
            }
        },
        "output" : {
            "body_file" : body_file,
            "control_point_file" : "dev/results/diamond_5_deg_full_{0}_control_points.vtk".format(grid)
        }
    }

    # Dump
    input_file = "dev/diamond_input.json"
    with open(input_file, 'w') as input_handle:
        json.dump(input_dict, input_handle, indent=4)

    # Run
    sp.run(["./machline.exe", input_file])
    
    # Run shock-expansion comparison
    airfoil = fl.DiamondAirfoil(5.0, 1.0)
    airfoil.set_state(M, alpha, gamma, p_inf, T_inf, c_inf, rho, 1.0e-5, 800.0, 0.7)
    p2, p3, p4, p5 = airfoil.get_pressures()
    x = 0.5*gamma*M**2
    Cp2 = (p2/p_inf-1.0)/x
    Cp3 = (p3/p_inf-1.0)/x
    Cp4 = (p4/p_inf-1.0)/x
    Cp5 = (p5/p_inf-1.0)/x
    CL, CD = airfoil.get_inviscid_lift_and_drag_coefs()
    Cx = -CL*np.sin(np.radians(alpha)) + CD*np.cos(np.radians(alpha))
    Cz = CL*np.cos(np.radians(alpha)) + CD*np.sin(np.radians(alpha))

    # Print table
    print("{0:<20}{1:<20}".format("Surface", "Cp"))
    print("".join(["-"]*40))
    print("{0:<20}{1:<20}".format(2, Cp2))
    print("{0:<20}{1:<20}".format(3, Cp3))
    print("{0:<20}{1:<20}".format(4, Cp4))
    print("{0:<20}{1:<20}".format(5, Cp5))
    print()
    print("CL:", CL)
    print("CD:", CD)
    print("Cx:", Cx)
    print("Cz:", Cz)

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
    save_loc = 'dev/results/diamond_5_deg_full_{0}_data.csv'.format(grid)
    pvs.SaveData(save_loc, proxy=plot, PointDataArrays=data_to_process, FieldAssociation='Point Data', Precision=12)

    # Read in data
    data = np.genfromtxt(save_loc, delimiter=',', skip_header=1)
    
    # Plot data from MachLine
    plt.figure()
    plt.plot(data[:,2], data[:,0], 'ks', label='2nd')
    plt.plot(data[:,2], data[:,1], 'kv', label='Ise.')

    # Plot data from shock-expansion theory
    x = np.linspace(0.0, 1.0, 100)
    Cp_upper = np.ones_like(x)
    Cp_upper[:50] *= Cp3
    Cp_upper[50:] *= Cp5
    Cp_lower = np.ones_like(x)
    Cp_lower[:50] *= Cp2
    Cp_lower[50:] *= Cp4
    plt.plot(x, Cp_upper, 'k--', label='SE Upper')
    plt.plot(x, Cp_lower, 'k-.', label='SE Lower')

    # Format
    plt.xlabel('$x$')
    plt.ylabel('$C_p$')
    plt.legend()
    plt.show()