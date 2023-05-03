import os
import json

import numpy as np
import paraview.simple as pvs
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, write_input_file, cases, line_styles, quad_labels
from studies.paraview_functions import extract_all_data, get_data_column_from_array, extract_plane_slice


RERUN_MACHLINE = True
study_dir = "studies/supersonic_love_delta_wing/"
M = 1.62
semispan_locs = [22.5, 64.1]
b_mid = (0.463/2/.230)
convergence_plot_densities = ['coarse', 'medium', 'fine']


def plot_pressure_slices(pressure_rule, angles_of_attack, mesh_density):

    # Get label for pressure rule
    if pressure_rule == "isentropic":
        pressure_rule_label = "C_p_ise"
    elif pressure_rule == "slender-body":
        pressure_rule_label = "C_p_sln"
    elif pressure_rule == "linear":
        pressure_rule_label = "C_p_lin"
    elif pressure_rule == "second-order":
        pressure_rule_label = "C_p_2nd"
    
    # Loop over angles of attack
    for j, alpha in enumerate(angles_of_attack):

        # Iterate over percent semispan locations
        for k, semi in enumerate(semispan_locs):

            # Iterate over cases
            for l, case in enumerate(cases):

                # Get MachLine results
                body_file = study_dir + "results/delta_wing_{0}_deg_{1}{2}.vtk".format(alpha, mesh_density, quad_labels[l])
                origin = [0.0, 0.0, semi*b_mid*0.01]
                normal = [0.0, 0.0, 1.0]
                headers, data = extract_plane_slice(body_file, normal, origin, filter=False, which_data='cell')
                x_ML = get_data_column_from_array(headers, data, "centroid:0")
                C_P_ML = get_data_column_from_array(headers, data, pressure_rule_label)

                # pull in experimental data based on Angle of Attack
                if alpha == 0.0:
                    surf = ['']
                elif alpha == 8.5:
                    surf = ['_lower']
                else:
                    surf = ['_upper', '_lower']
            
                # Initialize list to store marker shapes
                fill = ['k', 'w']

                # Iterate over upper and lower surface results from experimental data and plot each surface as a different marker
                plt.figure()
                for i, surface in enumerate(surf):
                    surf_label = "Exp. " + surface[1:].title()

                    data_file_exp = study_dir + "experimental_data/delta_wing_exp_{0}_deg_aoa_{1}_semispan{2}.csv".format(alpha, semi, surface)
                    experimental_data = np.genfromtxt(data_file_exp, delimiter=",", dtype=float) # x in column 0, Cp in column 1 from Love, page 57

                    plt.plot(experimental_data[:,0], experimental_data[:,1], 'ks', label=surf_label, mfc=fill[i])

                # Plot theoretical data at 0 degrees alpha as presented by Love's work
                if alpha == 0.0:
            
                    theory_loc = study_dir + "experimental_data/delta_wing_exp_theory_0.0_deg_aoa_{0}_semispan.csv".format(semi)
                    theory_data = np.genfromtxt(theory_loc, delimiter=",", dtype=float)

                    plt.plot(theory_data[:,0], theory_data[:,1], 'k:', label="Love-Hayes")
            
                # Differentiate between the upper and lower surfaces of MachLine's results
                mid = round(len(x_ML)/2)

                # Plot upper surface
                plt.plot(x_ML[:mid], C_P_ML[:mid], 'k-', label='ML Upper')
            
                # Plot lower surface
                plt.plot(x_ML[mid:], C_P_ML[mid:], 'k--', label = 'ML Lower')

                # Format plot
                plt.xlabel("$\\frac{x}{c_r}$")
                plt.ylabel("$C_P$")
                plt.gca().invert_yaxis()
                plt.legend()

                # Save figure
                if not os.path.exists(study_dir + 'plots/'):
                    os.makedirs(study_dir + 'plots/')
                plot_loc = study_dir + 'plots/delta_wing_comparison_{0}_semispan_{1}_alpha_{2}_{3}.pdf'.format(semi, alpha, pressure_rule, case)
                plt.savefig(plot_loc)
                plt.close()


def plot_force_convergence_over_alpha(alpha_list):

    # Pull in experimental data
    CD_exp_loc = study_dir + 'experimental_data/delta_wing_exp_CD.csv'
    CL_exp_loc = study_dir + 'experimental_data/delta_wing_exp_CL.csv'
    CM_exp_loc = study_dir + 'experimental_data/delta_wing_exp_CM.csv'

    CD_exp = np.genfromtxt(CD_exp_loc, delimiter=',')
    CL_exp = np.genfromtxt(CL_exp_loc, delimiter=',')
    CM_exp = np.genfromtxt(CM_exp_loc, delimiter=',')

    # Pull in CPanel data to compare against
    CPanel_CD_loc = study_dir + 'experimental_data/Davis_CPanel_CD.csv'
    CPanel_CL_loc = study_dir + 'experimental_data/Davis_CPanel_CL.csv'

    CPanel_CD = np.genfromtxt(CPanel_CD_loc, delimiter=',')
    CPanel_CL = np.genfromtxt(CPanel_CL_loc, delimiter=',')

    # Initialize lists to store force values over alpha range
    CD_coarse = []
    CL_coarse = []
    CD_medium = []
    CL_medium = []
    CD_fine = []
    CL_fine = []

    # Iterate over alpha list and mesh densities

    # Results for coarse mesh
    for alpha in alpha_list:
        results_loc = study_dir + 'results/delta_wing_{0}_coarse.json'.format(alpha)

        # Pull MachLine force data
        json_string = open(results_loc).read()
        json_vals = json.loads(json_string) 

        force_coefs = json_vals['total_forces']
        CD_coarse.append(force_coefs['Cx'] * np.cos(alpha*np.pi/180) + force_coefs['Cy'] * np.sin(alpha*np.pi/180))
        CL_coarse.append(-force_coefs['Cx'] * np.sin(alpha*np.pi/180) + force_coefs['Cy'] * np.cos(alpha*np.pi/180))

    # Results for medium mesh
    for alpha in alpha_list:
        results_loc = study_dir + 'results/delta_wing_{0}_medium.json'.format(alpha)

        # Pull MachLine force data
        json_string = open(results_loc).read()
        json_vals = json.loads(json_string) 

        force_coefs = json_vals['total_forces']
        CD_medium.append(force_coefs['Cx'] * np.cos(alpha*np.pi/180) + force_coefs['Cy'] * np.sin(alpha*np.pi/180))
        CL_medium.append(-force_coefs['Cx'] * np.sin(alpha*np.pi/180) + force_coefs['Cy'] * np.cos(alpha*np.pi/180))


    # Results for fine mesh
    for alpha in alpha_list:
        results_loc = study_dir + 'results/delta_wing_{0}_fine.json'.format(alpha)

        # Pull MachLine force data
        json_string = open(results_loc).read()
        json_vals = json.loads(json_string) 

        force_coefs = json_vals['total_forces']
        CD_fine.append(force_coefs['Cx'] * np.cos(alpha*np.pi/180) + force_coefs['Cy'] * np.sin(alpha*np.pi/180))
        CL_fine.append(-force_coefs['Cx'] * np.sin(alpha*np.pi/180) + force_coefs['Cy'] * np.cos(alpha*np.pi/180))
    
    
    # Plot results against experimental data

    # CD plot
    plt.figure()
    plt.scatter(alpha_list, CD_coarse, label="MachLine: coarse", marker='o', s= 300, edgecolors='k', facecolors='none')
    plt.scatter(alpha_list, CD_medium, label="MachLine: medium", marker='o', s= 100, edgecolors='k', facecolors='none')
    plt.scatter(alpha_list, CD_fine, label="MachLine: fine", marker='o', s= 30, edgecolors='k', facecolors='none')
    plt.scatter(CD_exp[:,0], CD_exp[:,1], label='Experimental Data', marker='s', edgecolors='k', facecolors='none')
    plt.scatter(CPanel_CD[:,0], CPanel_CD[:,1], label='CPanel', marker='D', edgecolors='k', facecolors='none')

    plt.legend()
    plt.xlabel('Angle of Attack [deg]')
    plt.ylabel('CD')
    plot_loc = study_dir + 'plots/delta_wing_CD_convergence.pdf'
    plt.savefig(plot_loc)
    plt.close()

    # CL plot
    plt.figure()
    plt.scatter(alpha_list, CL_coarse, label="MachLine: coarse", marker='o', s=300, edgecolors='k', facecolors='none')
    plt.scatter(alpha_list, CL_medium, label="MachLine: medium", marker='o', s=100, edgecolors='k', facecolors='none')
    plt.scatter(alpha_list, CL_fine, label="MachLine: fine", marker='o', s=30, edgecolors='k', facecolors='none')
    plt.scatter(CL_exp[:,0], CL_exp[:,1], label='Experimental Data', marker='s', edgecolors='k', facecolors='none')
    plt.scatter(CPanel_CL[:,0], CPanel_CL[:,1], label='CPanel', marker='D', edgecolors='k', facecolors='none')
    plt.legend()
    plt.xlabel('Angle of Attack [deg]')
    plt.ylabel('CL')
    plot_loc = study_dir + 'plots/delta_wing_CL_convergence.pdf'
    plt.savefig(plot_loc)
    plt.close()

    # Print CD and CL results along with exerimental data, allowing for quanitification
    # of difference between MachLine and experimental resutls for paper purposes.

    csv_loc = study_dir + 'error_quantification_test.csv'
    # Transpose all arrays to format for column vectors
    angles_col = np.array([alpha_list]).T
    CD_col = np.array([CD_fine]).T
    CD_ML = np.array([CD_fine]).T
    CL_ML = np.array([CL_fine]).T
    
    force_array = np.hstack((angles_col, CD_ML))
    force_array = np.hstack((force_array, CL_ML))
    # file = np.savetxt(csv_loc, force_array, delimiter=',')


def plot_force_alpha(alpha_list, M):
    # uses the fine mesh

    # Pull in experimental data
    CD_exp_loc = study_dir + 'experimental_data/delta_wing_exp_CD.csv'
    CL_exp_loc = study_dir + 'experimental_data/delta_wing_exp_CL.csv'
    CM_exp_loc = study_dir + 'experimental_data/delta_wing_exp_CM.csv'

    CD_exp = np.genfromtxt(CD_exp_loc, delimiter=',')
    CL_exp = np.genfromtxt(CL_exp_loc, delimiter=',')
    CM_exp = np.genfromtxt(CM_exp_loc, delimiter=',')

    # Pull in CPanel data to compare against
    CPanel_CD_loc = study_dir + 'experimental_data/Davis_CPanel_CD.csv'
    CPanel_CL_loc = study_dir + 'experimental_data/Davis_CPanel_CL.csv'
    CPanel_CD = np.genfromtxt(CPanel_CD_loc, delimiter=',')
    CPanel_CL = np.genfromtxt(CPanel_CL_loc, delimiter=',')

    # Initialize lists to store force values over alpha range
    CD_coarse = []
    CL_coarse = []
    CD_medium = []
    CL_medium = []
    CD_fine = []
    CL_fine = []

    # Results for fine mesh
    for alpha in alpha_list:

        # Pull MachLine force data
        results_loc = study_dir + 'results/delta_wing_{0}_fine.json'.format(alpha)
        with open(results_loc) as report_handle:
            report_dict = json.load(report_handle)

        force_coefs = report_dict['total_forces']
        CD_fine.append(force_coefs['Cx'] * np.cos(alpha*np.pi/180) + force_coefs['Cy'] * np.sin(alpha*np.pi/180))
        CL_fine.append(-force_coefs['Cx'] * np.sin(alpha*np.pi/180) + force_coefs['Cy'] * np.cos(alpha*np.pi/180))

        # Get which pressure rule was used for the forces
        pressure_for_forces = report_dict["input"]["post_processing"].get('pressure_for_forces', 'isentropic')
    
    # Plot results against experimental data

    # CD plot
    plt.figure()
    plt.plot(alpha_list, CD_fine, 'ko', label="MachLine", markersize=3)
    plt.plot(CD_exp[:,0], CD_exp[:,1], 'ks', label='Experiment', markersize=3)
    plt.plot(CPanel_CD[:,0], CPanel_CD[:,1], 'kv', label='CPanel', markersize=3)
    plt.legend()
    plt.xlabel('$\\alpha [^\circ]$')
    plt.ylabel('$C_D$')
    plt.ylim(bottom=0.0)
    plot_loc = study_dir + 'plots/delta_wing_CD_comparison_{0}.pdf'.format(pressure_for_forces)
    plt.savefig(plot_loc)
    plt.close()

    # CL plot
    plt.figure()
    plt.plot(alpha_list, CL_fine, 'ko', label="MachLine", markersize=3)
    plt.plot(CL_exp[:,0], CL_exp[:,1], 'ks', label='Experiment', markersize=3)
    plt.plot(CPanel_CL[:,0], CPanel_CL[:,1], 'kv', label='CPanel', markersize=3)
    plt.legend()
    plt.xlabel('$\\alpha [^\circ]$')
    plt.ylabel('$C_L$')
    plot_loc = study_dir + 'plots/delta_wing_CL_comparison_{0}.pdf'.format(pressure_for_forces)
    plt.savefig(plot_loc)
    plt.close()

    # Print CD and CL results along with exerimental data, allowing for quanitification
    # of difference between MachLine and experimental results for paper purposes.
    csv_loc = study_dir + 'error_quantification_test.csv'

    # Transpose all arrays to format for column vectors
    angles_col = np.array([alpha_list]).T
    CD_col = np.array([CD_fine]).T
    CD_ML = np.array([CD_fine]).T
    CL_ML = np.array([CL_fine]).T
    
    force_array = np.hstack((angles_col, CD_ML))
    force_array = np.hstack((force_array, CL_ML))
    # file = np.savetxt(csv_loc, force_array, delimiter=',')


def run_pressure_distribution_comparison():

    # Parameters
    angles_of_attack = [0.0, 2.0, 4.1, 8.5, 10.75]
    mesh_density = 'coarse' # coarse, medium, fine are the options created with 20, 40, and 80 nodes respectively

    # Run MachLine
    run_machline_cases(angles_of_attack, M, mesh_density)

    #for alpha in angles_of_attack:

    #    # Read into ParaView
    #    body_file = study_dir + "results/delta_wing_{0}_deg_{1}.vtk".format(alpha, mesh_density)
    #    data_reader = pvs.LegacyVTKReader(registrationName=body_file.replace(study_dir + "results", ""), FileNames=body_file)

    #    # Filter cell data to point data
    #    filter = pvs.CellDatatoPointData(registrationName='Filter', Input=data_reader)
    #    data_to_process = ['C_p_ise', 'C_p_2nd', 'C_p_lin', 'C_p_sln', 'centroid']
    #    filter.CellDataArraytoprocess = data_to_process
    #    fields_to_process = 'Point Data'

    #    # Iterate over semispan locations
    #    for i, percent_semispan in enumerate(semispan_locs):
    #        semispan_str = str(percent_semispan)
    #        percent_semi = percent_semispan * b_mid / 100.0

    #        # Slice mesh at each semispan location
    #        slicer = pvs.PlotOnIntersectionCurves(registrationName="Slice", Input=filter)
    #        slicer.SliceType = 'Plane'
    #        

    #        slicer.SliceType.Origin = origin
    #        slicer.SliceType.Normal = normal

    #        # Extract and save data
    #        plot = pvs.PlotData(registrationName="Plot")
    #        view = pvs.CreateView('XYChartView')
    #        display = pvs.Show(plot, view, 'XYChartRepresentation')
    #        display.XArrayName = 'centroid_X'
    #        view.Update()
    #        save_loc = study_dir + 'results/delta_wing_{0}_semispan_{1}_deg_results.csv'.format(semispan_str, alpha)
    #        pvs.SaveData(save_loc, proxy=plot, CellDataArrays=data_to_process, FieldAssociation='Point Data')

    # Plot pressure rule method over a range of angles of attack at each semispan location
    pressure_rule = "isentropic" # isentropic, second-order, slender-body, or linear
    plot_pressure_slices(pressure_rule, angles_of_attack, mesh_density)


def run_force_comparison():

    # Parameters
    angles_of_attack = [-6., -5., -4., -3., -2., -1., 0., 1., 2., 3., 4., 5.]
    mesh_density = 'fine' # coarse, medium, fine are the options created with 20, 40, and 80 nodes respectively

    # Run MachLine
    run_machline_cases(angles_of_attack, M, mesh_density)

    # Plot
    plot_force_alpha(angles_of_attack, M)
    plot_force_convergence_over_alpha(angles_of_attack, M)


def run_machline_cases(angles_of_attack, M, mesh_density):

    # Iterate over angles of attack
    for alpha in angles_of_attack:
        
        # Declare MachLine input
        body_file = study_dir + "results/delta_wing_{0}_deg_{1}.vtk".format(alpha, mesh_density)
        input_dict = {
            "flow": {
                "freestream_velocity": [np.cos(np.radians(alpha)), np.sin(np.radians(alpha)), 0.0],
                "gamma" : 1.4,
                "freestream_mach_number" : M
            },
            "geometry": {
                "file": study_dir + "meshes/delta_wing_clustered_mesh_{0}.vtk".format(mesh_density),
                "spanwise_axis" : "+z",
                "max_continuity_angle" : 1.0,
                "wake_model": {
                    "wake_present" : True,
                    "append_wake" : False
                },
                "reference": {
                    "area": b_mid*1.0 # 1 represents the root chord for this wing. b_mid is the semispan
                }
            },
            "solver": {
                "formulation": "morino",
                "matrix_solver": "GMRES",
                "run_checks": False
            },
            "post_processing" : {
                "pressure_rules" : {
                    "second-order" : True,
                    "isentropic" : True,
                    "slender-body" : True,
                    "linear" : True
                },
                "pressure_for_forces" : 'second-order'
            },
            "output" : {
                "body_file" :          body_file,
                "control_point_file" : study_dir + "results/delta_wing_{0}_deg_{1}_control_points.vtk".format(alpha,mesh_density),
                "wake_file": study_dir + "results/delta_wing_{0}_deg_{1}_wake.vtk".format(alpha, mesh_density),
                "report_file": study_dir + "reports/delta_wing_{0}_{1}.json".format(alpha, mesh_density)
            }
        }

        # Dump
        input_file = study_dir + "input_files/delta_wing_input.json"
        write_input_file(input_dict, input_file)

        # Run
        run_quad(input_file, run=RERUN_MACHLINE)
        

if __name__=="__main__":

    # Run comparisons
    run_pressure_distribution_comparison()
    #run_force_comparison()