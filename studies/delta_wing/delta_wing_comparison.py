import json
import subprocess as sp
import numpy as np
import paraview.simple as pvs
import matplotlib.pyplot as plt
from sys import exit
from os import listdir, chdir, getcwd

def data_plot(comp_method, angles_of_attack, semispan_locations):

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


            # Initialize figure and shapes to be plotted
            plt.figure(figsize=(6,4))
            shape = ["o", "^", "s", "v"]

            data_location = 'studies/delta_wing/results/delta_wing_{0}_semispan_{1}_deg_results.csv'.format(semi, AoA)

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

            # pull in experimental data based on Angle of Attack
            if AoA == 0.0:
                surf = ['']
            elif AoA == 8.5:
                surf = ['_lower']
            else:
                surf = ['_upper', '_lower']
            
            # Initialize list to store marker shapes
            shape = ['o', 's']

            
            # Iterate over upper and lower surface results from experimental data and plot each surface as a different marker
            for i,surface in enumerate(surf):
                surf_label = "Experimental " + surface[1:]
                
                exper_data_loc = "studies/delta_wing/experimental_data/delta_wing_exp_{0}_deg_aoa_{1}_semispan{2}.csv".format(AoA,semi,surface)
                experimental_data = np.genfromtxt(exper_data_loc, delimiter=",", dtype=float) # x in column 0, Cp in column 1 from Love, page 57
                
                plt.plot(experimental_data[:,0], experimental_data[:,1], color='k', label=surf_label, marker = shape[i], linestyle="none", fillstyle='full')
                
            # Plot theoretical data at 0 degrees AoA as presented by Love's work
            if AoA == 0.0:
            
                theory_loc = "studies/delta_wing/experimental_data/delta_wing_exp_theory_0.0_deg_aoa_{0}_semispan.csv".format(semi)
                theory_data = np.genfromtxt(theory_loc, delimiter=",", dtype=float)

                plt.plot(theory_data[:,0], theory_data[:,1], label="Theory", linewidth=1.5,color='k')
            
            # Differentiate between the upper and lower surfaces of MachLine's results
            half = round(len(data[:,x_loc])/2)

            # breakpoint()
            # Plot upper surface
            plt.plot(data[0:half,x_loc]/max(data[:,x_loc]), data[:half,Cp_loc],color='k',linestyle='--',label = 'MachLine upper',marker = '.',markersize=5,fillstyle='full')
            
            # Plot lower surface
            plt.plot(data[half::,x_loc]/max(data[:,x_loc]),data[half:,Cp_loc],color='k',linestyle='--',label = 'MachLine lower',marker = '^', markersize=4,fillstyle='full')
            # plt.plot(data[:,x_loc]/max(data[:,x_loc]), data[:,Cp_loc],color='k', label = 'MachLine', marker = "s", linestyle="none")

            # Format plot
            plt.xlabel(r"$\frac{x}{c_{root}}$")
            y_title = r"$C_p$"
            plt.ylabel(y_title)
            plt.gca().invert_yaxis()
            plt.legend()
            # plt.title(f"{semi} percent semispan at {aoa_formatted}")

            # Save figure
            plot_loc = 'studies/delta_wing/plots/delta_wing_comparison_{0}_semispan_{1}_aoa.pdf'.format(semi, AoA)
            plt.savefig(plot_loc)
            plt.show()

def force_plot(AoA_list,M,densities):

    # Pull in experimental data
    CD_exp_loc = 'studies/delta_wing/experimental_data/delta_wing_exp_CD.csv'
    CL_exp_loc = 'studies/delta_wing/experimental_data/delta_wing_exp_CL.csv'
    CM_exp_loc = 'studies/delta_wing/experimental_data/delta_wing_exp_CM.csv'

    CD_exp = np.genfromtxt(CD_exp_loc, delimiter=',')
    CL_exp = np.genfromtxt(CL_exp_loc, delimiter=',')
    CM_exp = np.genfromtxt(CM_exp_loc, delimiter=',')

    # Pull in CPanel data to compare against
    CPanel_CD_loc = 'studies/delta_wing/experimental_data/Davis_CPanel_CD.csv'
    CPanel_CL_loc = 'studies/delta_wing/experimental_data/Davis_CPanel_CL.csv'

    CPanel_CD = np.genfromtxt(CPanel_CD_loc, delimiter=',')
    CPanel_CL = np.genfromtxt(CPanel_CL_loc, delimiter=',')

    # Initialize lists to store force values over AoA range
    CD_coarse = []
    CL_coarse = []
    CD_semi_fine = []
    CL_semi_fine = []
    CD_fine = []
    CL_fine = []

    # Iterate over AoA list and mesh densities

    # Results for coarse mesh
    for AoA in AoA_list:
        results_loc = 'studies/delta_wing/results/delta_wing_{0}_coarse.json'.format(AoA)

        # Pull MachLine force data
        json_string = open(results_loc).read()
        json_vals = json.loads(json_string) 

        loc = json_vals['total_forces']
        CD_coarse.append(loc['Cx'] * np.cos(AoA*np.pi/180) + loc['Cy'] * np.sin(AoA*np.pi/180))
        CL_coarse.append(-loc['Cx'] * np.sin(AoA*np.pi/180) + loc['Cy'] * np.cos(AoA*np.pi/180))

    # Results for semi_fine mesh
    for AoA in AoA_list:
        results_loc = 'studies/delta_wing/results/delta_wing_{0}_semi_fine.json'.format(AoA)

        # Pull MachLine force data
        json_string = open(results_loc).read()
        json_vals = json.loads(json_string) 

        loc = json_vals['total_forces']
        CD_semi_fine.append(loc['Cx'] * np.cos(AoA*np.pi/180) + loc['Cy'] * np.sin(AoA*np.pi/180))
        CL_semi_fine.append(-loc['Cx'] * np.sin(AoA*np.pi/180) + loc['Cy'] * np.cos(AoA*np.pi/180))


    # Results for fine mesh
    for AoA in AoA_list:
        results_loc = 'studies/delta_wing/results/delta_wing_{0}_fine.json'.format(AoA)

        # Pull MachLine force data
        json_string = open(results_loc).read()
        json_vals = json.loads(json_string) 

        loc = json_vals['total_forces']
        CD_fine.append(loc['Cx'] * np.cos(AoA*np.pi/180) + loc['Cy'] * np.sin(AoA*np.pi/180))
        CL_fine.append(-loc['Cx'] * np.sin(AoA*np.pi/180) + loc['Cy'] * np.cos(AoA*np.pi/180))
    
    
    # Plot results against experimental data

    # CD plot
    plt.figure()
    plt.scatter(AoA_list, CD_coarse, label="MachLine: coarse", marker='o', s= 300, edgecolors='k', facecolors='none')
    plt.scatter(AoA_list, CD_semi_fine, label="MachLine: medium", marker='o', s= 100, edgecolors='k', facecolors='none')
    plt.scatter(AoA_list, CD_fine, label="MachLine: fine", marker='o', s= 30, edgecolors='k', facecolors='none')
    plt.scatter(CD_exp[:,0], CD_exp[:,1], label='Experimental Data', marker='s', edgecolors='k', facecolors='none') # errorbars are not included due to webplotdigitizer uncertainty magnitude
    plt.scatter(CPanel_CD[:,0], CPanel_CD[:,1], label='CPanel', marker='D', edgecolors='k', facecolors='none')

    plt.legend()
    plt.xlabel('Angle of Attack [deg]')
    plt.ylabel('CD')
    plot_loc = 'studies/delta_wing/plots/delta_wing_CD_comparison.pdf'
    plt.savefig(plot_loc)
    plt.show()

    # CL plot
    plt.figure()
    plt.scatter(AoA_list, CL_coarse, label="MachLine: coarse", marker='o', s=300, edgecolors='k', facecolors='none')
    plt.scatter(AoA_list, CL_semi_fine, label="MachLine: medium", marker='o', s=100, edgecolors='k', facecolors='none')
    plt.scatter(AoA_list, CL_fine, label="MachLine: fine", marker='o', s=30, edgecolors='k', facecolors='none')
    plt.scatter(CL_exp[:,0], CL_exp[:,1], label='Experimental Data', marker='s', edgecolors='k', facecolors='none') # errorbars are not included due to webplotdigitizer uncertainty magnitude
    plt.scatter(CPanel_CL[:,0], CPanel_CL[:,1], label='CPanel', marker='D', edgecolors='k', facecolors='none')
    plt.legend()
    plt.xlabel('Angle of Attack [deg]')
    plt.ylabel('CL')
    plot_loc = 'studies/delta_wing/plots/delta_wing_CL_comparison.pdf'
    plt.savefig(plot_loc)
    plt.show()


    # Print CD and CL results along with exerimental data, allowing for quanitification
    # of difference between MachLine and experimental resutls for paper purposes.

    csv_loc = 'studies/delta_wing/error_quantification_test.csv'
    # Transpose all arrays to format for column vectors
    angles_col = np.array([AoA_list]).T
    CD_col = np.array([CD_fine]).T
    CD_ML = np.array([CD_fine]).T
    CL_ML = np.array([CL_fine]).T
    
    force_array = np.hstack((angles_col, CD_ML))
    force_array = np.hstack((force_array, CL_ML))
    # file = np.savetxt(csv_loc, force_array, delimiter=',')


if __name__=="__main__":

    # Parameters
    M = 1.62
    R_G = 287.058
    gamma = 1.4
    T_inf = 300.0
    c_inf = np.sqrt(gamma*R_G*T_inf)
    # angles_of_attack = [0.0, 2.0, 4.1, 8.5, 10.75]
    angles_of_attack = [-6., -5., -4., -3., -2., -1., 0., 1., 2., 3., 4., 5.]
    # angles_of_attack = [0.0]
    semispan_loc = [22.5, 64.1]
    b_half = (0.463/2/.230) # semispan length nondimensionalized by root chord
    R_A = 4.023
    # b_half = 0.2315 # semispan length for OpenVSP model
    # mesh_type = "clustered"
    mesh_type = 'clustered'
    mesh_density = 'semi_fine' # coarse, semi_fine, fine are the options created with 20, 40, and 80 nodes respectively
    mesh_filetype = 'vtk'
    force_plot_densities = ['coarse', 'semi_fine', 'fine']

    # Check working directory and re-route if necessary
    check_dir = getcwd()

    if 'delta_wing' in check_dir:
        chdir("../../")

    # Iterate over angles of attack
    for alpha in angles_of_attack:
        if 'VSP' in mesh_type:
            freestream = [M*c_inf*np.cos(np.radians(alpha)), 0.0, M*c_inf*np.sin(np.radians(alpha))]
        else:
            freestream = [M*c_inf*np.cos(np.radians(alpha)), M*c_inf*np.sin(np.radians(alpha)), 0.0]
        
        # Declare MachLine input
        body_file = "studies/delta_wing/results/delta_wing_{0}_deg_{1}.vtk".format(alpha,mesh_density)
        input_dict = {
            "flow": {
                "freestream_velocity": freestream,
                "gamma" : gamma,
                "freestream_mach_number" : M
            },
            "geometry": {
                "file": "studies/delta_wing/meshes/delta_wing_{0}_mesh_{1}.{2}".format(mesh_type, mesh_density, mesh_filetype),
                "wake_model": {
                    "wake_present" : True,
                    "append_wake" : True,
                    "trefftz_distance": 20,
                },
                "reference": {
                    "area": (b_half) * 1 # 1 represents the root chord for this wing. b_half is the semispan
                }
            },
            "solver": {
                "formulation": "morino",
                # "matrix_solver": "QRUP",
                "run_checks": True
            },
            "post_processing" : {
                "pressure_rules" : {
                    "second-order" : True,
                    "isentropic" : True,
                    "slender-body" : True,
                    "linear" : True
                },
                # "pressure_for_forces" : 'slender-body'
                
            },
            "output" : {
            "verbose": True,
            "body_file" :          body_file,
            "control_point_file" : "studies/delta_wing/results/delta_wing_{0}_deg_{1}_control_points.vtk".format(alpha,mesh_density),
            "wake_file": "studies/delta_wing/results/delta_wing_{0}_deg_{1}_wake.vtk".format(alpha, mesh_density),
            "report_file": "studies/delta_wing/results/delta_wing_{0}_{1}.json".format(alpha, mesh_density)
            }
        }

        # Dump
        input_file = "studies/delta_wing/input_files/delta_wing_input.json"
        with open(input_file, 'w') as input_handle:
            json.dump(input_dict, input_handle, indent=4)

        # Run
        # sp.run(["./machline.exe", input_file])
        
        # Verify that MachLine execution was successful
        control_point_file = "delta_wing_{0}_deg_{1}_control_points.vtk".format(alpha,mesh_density)
        if control_point_file not in listdir("studies/delta_wing/results/"):
            print(f"The desired file cannot be located:  studies/delta_wing/results/{control_point_file}")
            exit()

        
        # Read into ParaView
        data_reader = pvs.LegacyVTKReader(registrationName=body_file.replace("studies/delta_wing/results", ""), FileNames=body_file)

        # Filter cell data to point data
        filter = pvs.CellDatatoPointData(registrationName='Filter', Input=data_reader)
        data_to_process = ['C_p_ise', 'C_p_2nd', 'C_p_lin', 'C_p_sln', 'centroid']
        filter.CellDataArraytoprocess = data_to_process
        fields_to_process = 'Point Data'

        # Iterate over semispan locations
        for i, percent_semispan in enumerate(semispan_loc):
            semispan_str = str(percent_semispan)
            percent_semi = percent_semispan * b_half / 100.0

            # Slice mesh at each semispan location
            slicer = pvs.PlotOnIntersectionCurves(registrationName= "Slice", Input=filter)
            slicer.SliceType = 'Plane'
            # slicer.HyperTreeGridSlicer = 'Plane'
            
            if 'VSP' in mesh_type:
                origin = [0.0, percent_semi, 0.0]
                normal = [0.0, 1.0, 0.0]
            else:
                origin = [0.0, 0.0, percent_semi]
                normal = [0.0, 0.0, 1.0]

            slicer.SliceType.Origin = origin
            slicer.SliceType.Normal = normal


            # Extract and save data
            plot = pvs.PlotData(registrationName="Plot")
            view = pvs.CreateView('XYChartView')
            display = pvs.Show(plot, view, 'XYChartRepresentation')
            display.XArrayName = 'centroid_X'
            view.Update()
            save_loc = 'studies/delta_wing/results/delta_wing_{0}_semispan_{1}_deg_results.csv'.format(semispan_str, alpha)
            pvs.SaveData(save_loc, proxy=plot, ChooseArraysToWrite=1, CellDataArrays=data_to_process, FieldAssociation=fields_to_process)
    # Plot single computational method over a range of angles of attack at each semispan location
    computational_method = "isentropic" # isentropic, second order, slender-body, or linear
    # data_plot(computational_method, angles_of_attack, semispan_loc)

    # Plot CL CD CM
    force_plot(angles_of_attack, M, force_plot_densities)