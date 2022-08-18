import json
import subprocess as sp
import numpy as np
import paraview.simple as pvs
import matplotlib.pyplot as plt
from sys import exit
from os import listdir, chdir, getcwd

def data_plot(comp_method, angles_of_attack, semispan_locations, semispan_length):

    # Reformat computational method
    if comp_method == "isentropic":
        comp_type = "C_p_ise"
    elif comp_method == "slender-body":
        comp_type = "C_p_sln"
    elif comp_method == "linear":
        comp_type = "C_p_lin"
    elif comp_method == "second order":
        comp_type = "C_p_2nd"

    
    # Loop over angles of attack
    for j,AoA in enumerate(angles_of_attack):

        # Iterate over percent semispan locations
        for k, semi in enumerate(semispan_locations):

            # pull in experimental data
            exper_data_loc = "studies/delta_wing/experimental_data/delta_wing_exp_{0}_deg_aoa_{1}_semispan.csv".format(AoA,semi)
            experimental_data = np.genfromtxt(exper_data_loc, delimiter=",", dtype=float) # x in column 0, Cp in column 1 from Love, page 57

            # Initialize figure and shapes to be plotted
            plt.figure(figsize=(10,8))
            shape = ["o", "^", "s", "v"]

            data_location = 'studies/delta_wing/delta_wing_{0}_semispan_{1}_deg_results.csv'.format(semi, AoA)

            # Identify which column of data is associated with selected comp method
            column_data = np.genfromtxt(data_location, delimiter=",", dtype=str)
            column_data = np.char.replace(column_data[0,:],'"',"")

            # Identify locations of method data selected by user
            Cp_loc = -1
            x_loc = -1
            for i,col in enumerate(column_data):
                if col == comp_type:
                    Cp_loc = i
                if col == "centroid:0":
                    x_loc = i
                # Break out of loop if all points have been located
                if min(Cp_loc,x_loc) > -1:
                    break
            
            # Read in data without column headers
            data = np.genfromtxt(data_location, delimiter=",", skip_header=1, dtype=float)

            # Format angle of attack for legend
            aoa_formatted = str(AoA) + r"$^{\circ}$ deg"

            # Plot data
            plt.plot(data[:,x_loc]/max(data[:,x_loc]), data[:,Cp_loc],color='k', label = 'MachLine')
            # plt.plot(data[:,x_loc]/max(data[:,x_loc]), data[:,Cp_loc],color='k', label = 'MachLine', marker = "s", linestyle="none")
            # plt.plot(experimental_data[:,0], experimental_data[:,1], color='k', marker=".", label="Experimental", linestyle="none", fillstyle="full")
            plt.plot(experimental_data[:,0], experimental_data[:,1], color='k', label="Experimental", marker = "o", linestyle="none", fillstyle='full')


            # Format plot
            plt.xlabel(r"$\frac{x}{l}$")
            y_title = r"$C_p$ " + comp_method
            plt.ylabel(y_title)
            plt.gca().invert_yaxis()
            plt.legend()
            plt.title(f"{semi} percent semispan at {aoa_formatted}")

            # Save figure
            plot_loc = 'studies/delta_wing/plots/delta_wing_comparison_{0}_semispan_{1}_aoa.pdf'.format(semi, AoA)
            plt.savefig(plot_loc)
            plt.show()

if __name__=="__main__":

    # Parameters
    M = 1.62
    R_G = 287.058
    gamma = 1.4
    T_inf = 300.0
    c_inf = np.sqrt(gamma*R_G*T_inf)
    angles_of_attack = [0.0]#, 2.0, 4.1, 8.5, 10.75]
    # angles_of_attack = [0.0]
    semispan_loc = [22.5, 64.1]
    b_half = .2315 # meters
    mesh_density = "coarse_test"

    # Check working directory and re-route if necessary
    check_dir = getcwd()

    if 'delta_wing' in check_dir:
        chdir("../../")

    # Iterate over angles of attack
    for alpha in angles_of_attack:
        
        # Declare MachLine input
        body_file = "studies/delta_wing/delta_wing_{0}_deg_{1}.vtk".format(alpha,mesh_density)
        input_dict = {
            "flow": {
                "freestream_velocity": [M*c_inf*np.cos(np.radians(alpha)), 0.0, M*c_inf*np.sin(np.radians(alpha))],
                "gamma" : gamma,
                "freestream_mach_number" : M
            },
            "geometry": {
                "file": "studies/delta_wing/meshes/delta_wing_{0}_mesh.vtk".format(mesh_density),
                "spanwise_axis" : "+y",
                "mirror_about": "xy",
                "wake_model": {
                    "wake_present" : True,
                    "append_wake" : False
                },
                "reference": {
                    "area": 1.0
                }
            },
            "solver": {
                "formulation": "morino",
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
            "body_file" :          body_file,
            "mirrored_body_file":   "studies/delta_wing/delta_wing_{0}_deg_{1}_mirrored_body.vtk".format(alpha,mesh_density),
            "control_point_file" : "studies/delta_wing/delta_wing_{0}_deg_{1}_control_points.vtk".format(alpha,mesh_density),
            "mirrored_control_point_file": "studies/delta_wing/delta_wing_{0}_deg_{1}_mirrored_control_points.vtk".format(alpha,mesh_density)
            }
        }

        # Dump
        input_file = "studies/delta_wing/input_files/delta_wing_input.json"
        with open(input_file, 'w') as input_handle:
            json.dump(input_dict, input_handle, indent=4)

        # Run
        sp.run(["./machline.exe", input_file])
        
        # Verify that MachLine execution was successful
        control_point_file = "delta_wing_{0}_deg_{1}_control_points.vtk".format(alpha,mesh_density)
        if control_point_file not in listdir("studies/delta_wing/"):
            print(f"The desired file cannot be located:  studies/delta_wing/{control_point_file}")
            exit()



        # Read into ParaView
        data_reader = pvs.LegacyVTKReader(registrationName=body_file.replace("studies/", ""), FileNames=body_file)

        # Filter cell data to point data
        # filter = pvs.CellDatatoPointData(registrationName='Filter', Input=data_reader)
        data_to_process = ['C_p_ise', 'C_p_2nd', 'C_p_lin', 'C_p_sln', 'centroid']
        # filter.CellDataArraytoprocess = data_to_process

        # Iterate over semispan locations
        for i, percent_semispan in enumerate(semispan_loc):
            semispan_str = str(percent_semispan)
            percent_semispan = percent_semispan * b_half / 100.0

            # Slice mesh at each semispan location
            slicer = pvs.Slice(registrationName= "Slice", Input=data_reader)
            slicer.SliceType = 'Plane'
            # slicer.HyperTreeGridSlicer = 'Plane'
            slicer.SliceOffsetValues = [0.0]
            slicer.SliceType.Origin = [0.0, percent_semispan, 0.0]
            slicer.SliceType.Normal = [0.0, 1.0, 0.0]
            # slicer.HyperTreeGridSlicer = [0.0, 1.0, 0.0]

            # Extract and save data
            plot = pvs.PlotData(registrationName="Plot")
            view = pvs.CreateView('XYChartView')
            display = pvs.Show(plot, view, 'XYChartRepresentation')
            display.XArrayName = 'centroid_X'
            view.Update()
            save_loc = 'studies/delta_wing/delta_wing_{0}_semispan_{1}_deg_results.csv'.format(semispan_str, alpha)
            pvs.SaveData(save_loc, proxy=plot, ChooseArraysToWrite=1, CellDataArrays=data_to_process, FieldAssociation='Cell Data')

    # Plot single computational method over a range of angles of attack at each semispan location
    computational_method = "isentropic" # isentropic, second order, slender-body, or linear
    data_plot(computational_method, angles_of_attack, semispan_loc, b_half)