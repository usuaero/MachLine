import json
import subprocess as sp
import numpy as np
import paraview.simple as pvs
import matplotlib.pyplot as plt
from sys import exit
from os import listdir, chdir, getcwd

def comparison_plots(data_location, Mach):
    # Initialize location holders and other information for plotting
    sec_loc = -1
    ise_loc = -1
    lin_loc = -1
    sln_loc = -1
    x_loc = -1
    methods = []
    locations = []
    shape = [".", "v", "s", "^"]
    
    # Identify which column of data is associated with each calculation method
    column_data = np.genfromtxt(data_location, delimiter=",", dtype=str)
    column_data = np.char.replace(column_data[0,:],'"',"")

    # Identify locations of method data selected by user
    for i,col in enumerate(column_data):
        if col == "C_p_ise":
            ise_loc = i
        if col == "C_p_2nd":
            sec_loc = i
            methods.append("Second Order")
            locations.append(i)
        if col == "C_p_lin":
            lin_loc = i
            methods.append("Linear")
            locations.append(i)
        if col == "C_p_sln":
            sln_loc = i
            methods.append("Slender")
            locations.append(i)
        if col == "Points:0":
            x_loc = i

        # Break out of loop if all points have been located
        if min(sec_loc,ise_loc,lin_loc,sln_loc,x_loc) > -1:
            break

    # Verify that the isentropic rule was selected for comparison purposes
    if ise_loc == -1:
        print("In order to compare the computational methods, the isentropic rule must be selected. It will be used as the baseline method. Quitting...")
        exit()


    # Read in data without column headers
    data = np.genfromtxt(data_location, delimiter=",", skip_header=1, dtype=float)

    # Set baseline data that the other methods will be compared against
    baseline = data[:,ise_loc] # Isentropic data set as the baseline
    x_axis = data[:,x_loc] / max(data[:,x_loc]) # non-dimensionalized by total length

    # Plot data from MachLine
    plt.figure(figsize=(6,4))
    plt.yscale("log")

    # Differentiate positive and negative restults and plot for each identified method
    for i, method in enumerate(methods):
        # Initialize lists
        pos_diff = []
        pos_x_axis = []
        neg_diff = []
        neg_x_axis = []

        # Calculate method differences from the baseline
        diff = data[:,locations[i]] - baseline

        # Separate positive and negative values
        for j, val in enumerate(diff):
            if val >= 0.:
                pos_diff.append(val)
                pos_x_axis.append(x_axis[j])
            else:
                neg_diff.append(abs(val))
                neg_x_axis.append(x_axis[j])

        # Plot absolute value differences as different opacities
        if len(pos_diff) != 0:
            plt.plot(pos_x_axis, pos_diff, color='k', marker = shape[i], alpha=1, label=method, linestyle="none")
        if len(neg_diff) != 0:
            plt.plot(neg_x_axis, neg_diff, color='k', marker = shape[i], alpha=0.3, linestyle="none")


    # Format
    plt.xlabel(r'$\frac{x}{l}$')
    plt.ylabel('$Cp$ Difference')
    plt.legend()

    # Save figure
    plot_loc = 'studies/sears_haack/plots/sears_haack_M_{0}_comp_comparison.pdf'.format(Mach)
    plt.savefig(plot_loc)

    plt.show()

def data_plot(comp_method, Mach_numbers):

    # Reformat computational method
    if comp_method == "isentropic":
        comp_type = "C_p_ise"
    elif comp_method == "slender-body":
        comp_type = "C_p_sln"
    elif comp_method == "linear":
        comp_type = "C_p_lin"
    elif comp_method == "second order":
        comp_type = "C_p_2nd"

    # Initialize figure and shapes to be plotted
    shape = [".", "s", "^", "v"]

    # Loop over Mach Numbers
    for k,Mach in enumerate(Mach_numbers):
        plt.figure(figsize=(6,4))
        data_location = 'studies/sears_haack/machline_results/sears_haack_{0}_data_mach_{1}.csv'.format(grid,Mach)

        comparison_data_loc = 'studies/sears_haack/sears_haack_Stivers_mach_{}.csv'.format(Mach)
        
        # Identify which column of data is associated with selected comp method
        column_data = np.genfromtxt(data_location, delimiter=",", dtype=str)
        column_data = np.char.replace(column_data[0,:],'"',"")
        comparison_data = np.genfromtxt(comparison_data_loc, delimiter=',', dtype=float)

        # Identify locations of method data selected by user
        Cp_loc = -1
        x_loc = -1
        for i,col in enumerate(column_data):
            if col == comp_type:
                Cp_loc = i
            if col == "Points:0":
                x_loc = i
            # Break out of loop if all points have been located
            if min(Cp_loc,x_loc) > -1:
                break
        
        # Read in data without column headers
        data = np.genfromtxt(data_location, delimiter=",", skip_header=1, dtype=float)

        # Plot data
        plt.plot(data[:,x_loc]/max(data[:,x_loc]), data[:,Cp_loc],color='k', label = "MachLine", marker = '.', linestyle="none", fillstyle='full')
        plt.plot(comparison_data[:,0], comparison_data[:,1], 'k-', label="Stivers")

        # Format plot
        plt.xlabel("x/l")
        y_title = r"$C_p$"
        plt.ylabel(y_title)
        plt.gca().invert_yaxis()
        plt.legend()
        # plt.title(f"Mach {Mach}")

        # Save figure
        plot_loc = 'studies/sears_haack/plots/sears_haack_M_{0}.pdf'.format(Mach)
        plt.savefig(plot_loc)

        plt.show()

if __name__=="__main__":

    # Parameters
    M = 2
    Mach_Numbers = [2, 3] # Only for plotting purposes. The data needs to already be calculated previously
    R_G = 287.058
    gamma = 1.4
    T_inf = 300.0
    c_inf = np.sqrt(gamma*R_G*T_inf)
    alpha = 0.0
    grid = 9800 # 9800

    # Ensure correct working directory
    check_dir = getcwd()
    if 'sears_haack' in check_dir:
        chdir('../..')

    # Declare MachLine input
    body_file = "studies/sears_haack/machline_results/sears_haack_{0}_data_mach_{1}.vtk".format(grid, M)
    input_dict = {
        "flow": {
            "freestream_velocity": [M*c_inf*np.cos(np.radians(alpha)), 0.0, M*c_inf*np.sin(np.radians(alpha))],
            "gamma" : gamma,
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": "studies/sears_haack/meshes/sears_haack_mesh_love.tri",
            "spanwise_axis" : "+y",
            "wake_model": {
                "wake_present" : False,
            },
            "reference": {
                "area": 1.0
            }
        },
        "solver": {
            "formulation": "source-free",
            "run_checks": True,
            "control_point_offset": 1.1e-8
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
            "verbose": True,
            "body_file" : body_file,
            "control_point_file" : "studies/sears_haack/machline_results/sears_haack_{0}_data_mach_{1}_control_points.vtk".format(grid,M)
        }
    }

    # Dump
    input_file = "studies/sears_haack/sears_haack_input.json"
    with open(input_file, 'w') as input_handle:
        json.dump(input_dict, input_handle, indent=4)

    # Run
    sp.run(["./machline.exe", input_file])
    
    # Verify that MachLine execution was successful
    control_point_file = "sears_haack_{0}_data_mach_{1}_control_points.vtk".format(grid,M)
    if control_point_file not in listdir("studies/sears_haack/machline_results/"):
        print(f"The desired file cannot be located:  studies/sears_haack/{control_point_file}")
        exit()


    # Read into ParaView
    data_reader = pvs.LegacyVTKReader(registrationName=body_file.replace("studies/sears_haack/machline_results/", ""), FileNames=body_file)

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
    save_loc = 'studies/sears_haack/machline_results/sears_haack_{0}_data_mach_{1}.csv'.format(grid,M)
    pvs.SaveData(save_loc, proxy=plot, PointDataArrays=data_to_process, FieldAssociation='Point Data', Precision=12)
    # Read from json file which rules were selected
    json_string = open(input_file).read()
    json_vals = json.loads(json_string)

    # Initialize booleans
    second_order = json_vals["post_processing"]["pressure_rules"]["second-order"]
    isentropic = json_vals["post_processing"]["pressure_rules"]["isentropic"]
    slender = json_vals["post_processing"]["pressure_rules"]["slender-body"]
    linear = json_vals["post_processing"]["pressure_rules"]["linear"]


    # Call function to compare computational methods
    # comparison_plots(save_loc, M)

    # Plot single computational method over a range of mach numbers
    computational_method = "slender-body" # isentropic, second order, slender-body, or linear
    data_plot(computational_method,Mach_Numbers)