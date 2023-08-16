import os
import json

import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, write_input_file, cases, line_styles, quad_labels, get_order_of_convergence
from studies.paraview_functions import get_data_column_from_array, extract_plane_slice


RERUN_MACHLINE = False
PRESSURE_FOR_FORCES = 'isentropic'
study_dir = "studies/supersonic_love_delta_wing/"
results_dir = study_dir + "results/"
reports_dir = study_dir + "reports/"
plot_dir = study_dir + "plots/"
M = 1.62
semispan_locs = [22.5, 64.1]
b_mid = (0.463/2/.230)
mesh_densities = ['ultra_coarse', 'coarse', 'medium', 'fine', 'ultra_fine']
#mesh_densities = ["ultra_fine"]


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

                    plt.plot(theory_data[1:,0], theory_data[1:,1], 'k:', label="Love-Hayes")
            
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
                if case=='MH':
                    plt.legend()

                # Save figure
                if not os.path.exists(study_dir + 'plots/'):
                    os.makedirs(study_dir + 'plots/')
                plot_loc = study_dir + 'plots/delta_wing_comparison_{0}_{1}_semispan_{2}_alpha_{3}_{4}'.format(mesh_density, semi, alpha, pressure_rule, case)
                plt.savefig(plot_loc+'.pdf')
                plt.savefig(plot_loc+'.svg')
                plt.close()


def plot_force_convergence_over_alpha(alpha_list):

    # Pull in experimental data
    CD_exp = np.genfromtxt(study_dir + 'experimental_data/delta_wing_exp_CD.csv', delimiter=',')
    CL_exp = np.genfromtxt(study_dir + 'experimental_data/delta_wing_exp_CL.csv', delimiter=',')
    CM_exp = np.genfromtxt(study_dir + 'experimental_data/delta_wing_exp_CM.csv', delimiter=',')

    # Pull in CPanel data to compare against
    CPanel_CD = np.genfromtxt(study_dir + 'experimental_data/Davis_CPanel_CD.csv', delimiter=',')
    CPanel_CL = np.genfromtxt(study_dir + 'experimental_data/Davis_CPanel_CL.csv', delimiter=',')

    # Initialize storage
    CD = np.zeros((len(cases), len(mesh_densities), len(alpha_list)))
    CL = np.zeros((len(cases), len(mesh_densities), len(alpha_list)))
    Cx = np.zeros((len(cases), len(mesh_densities), len(alpha_list)))
    Cy = np.zeros((len(cases), len(mesh_densities), len(alpha_list)))
    l_avg = np.zeros(len(mesh_densities))

    # Iterate over cases
    for i, (case, quad_label) in enumerate(zip(cases, quad_labels)):

        # Iterate over mesh densities
        for j, mesh_density in enumerate(mesh_densities):

            # Iterate over alpha
            for k, alpha in enumerate(alpha_list):

                # Get trig ratios
                S = np.sin(np.radians(alpha))
                C = np.cos(np.radians(alpha))

                # Get report
                report_loc = reports_dir + 'delta_wing_{0}_{1}{2}.json'.format(alpha, mesh_density, quad_label)
                with open(report_loc, 'r') as results_handle:
                    report = json.load(results_handle) 

                # Store force data
                force_coefs = report['total_forces']
                Cx[i,j,k] = force_coefs['Cx']
                Cy[i,j,k] = force_coefs['Cy']
                CD[i,j,k] = force_coefs['Cx']*C + force_coefs['Cy']*S
                CL[i,j,k] = -force_coefs['Cx']*S + force_coefs['Cy']*C

                # Store average characteristic length
                l_avg[j] = report["mesh_info"]["average_characteristic_length"]
    
        # Plot convergence against experimental data

        # CD convergence plot
        err = np.abs((CD[i,:-1,:] - CD[i,-1,:])/CD[i,-1,:])
        plt.figure()
        for j, mesh_density in enumerate(mesh_densities[:-1]):
            plt.plot(alpha_list, err[j,:], 'o', label=mesh_density.replace("_", " ").title(), markersize=10-3*j, mec='k', mfc='none')
        #for j, mesh_density in enumerate(mesh_densities):
        #    plt.plot(alpha_list, CD[i,j,:], 'o', label=mesh_density.replace("_", " ").title(), markersize=10-2*j, mec='k', mfc='none')
        #if case=="MH":
        #    plt.legend()
        plt.xlabel('$\\alpha\,[^\circ]$')
        plt.ylabel('Fractional Error in $C_D$')
        plt.yscale('log')
        plt.ylim(bottom=1e-5, top=1e-2)
        plot_loc = study_dir + 'plots/delta_wing_CD_convergence{0}'.format(quad_label)
        plt.savefig(plot_loc + ".pdf")
        plt.savefig(plot_loc + ".svg")
        plt.close()

        # CL convergence plot
        err = np.abs((CL[i,:-1,:] - CL[i,-1,:])/CL[i,-1,:])
        ind = np.where(np.array(alpha_list) != 0.0)
        plt.figure()
        for j, mesh_density in enumerate(mesh_densities[:-1]):
            plt.plot(np.array(alpha_list)[ind], err[j,ind].flatten(), 'o', label=mesh_density.replace("_", " ").title(), markersize=10-3*j, mec='k', mfc='none')
        #for j, mesh_density in enumerate(mesh_densities):
        #    plt.plot(alpha_list, CL[i,j,:], 'o', label=mesh_density.replace("_", " ").title(), markersize=10-2*j, mec='k', mfc='none')
        #if case=="MH":
        #    plt.legend()
        plt.xlabel('$\\alpha\,[^\circ]$')
        plt.ylabel('Fractional Error in $C_L$')
        plt.yscale('log')
        plt.ylim(bottom=2e-4, top=1e-1)
        plot_loc = study_dir + 'plots/delta_wing_CL_convergence{0}'.format(quad_label)
        plt.savefig(plot_loc + ".pdf")
        plt.savefig(plot_loc + ".svg")
        plt.close()

        # Cx convergence plot
        err = np.abs((Cx[i,:-1,:] - Cx[i,-1,:])/Cx[i,-1,:])
        plt.figure()
        for j, mesh_density in enumerate(mesh_densities[:-1]):
            plt.plot(alpha_list, err[j,:], 'o', label=mesh_density.replace("_", " ").title(), markersize=10-3*j, mec='k', mfc='none')
        #for j, mesh_density in enumerate(mesh_densities):
        #    plt.plot(alpha_list, CD[i,j,:], 'o', label=mesh_density.replace("_", " ").title(), markersize=10-2*j, mec='k', mfc='none')
        #if case=="MH":
        #    plt.legend()
        plt.xlabel('$\\alpha\,[^\circ]$')
        plt.ylabel('Fractional Error in $C_x$')
        plt.yscale('log')
        plt.ylim(bottom=1e-5, top=1e-1)
        plot_loc = study_dir + 'plots/delta_wing_Cx_convergence{0}'.format(quad_label)
        plt.savefig(plot_loc + ".pdf")
        plt.savefig(plot_loc + ".svg")
        plt.close()

        # Cy convergence plot
        err = np.abs((Cy[i,:-1,:] - Cy[i,-1,:])/Cy[i,-1,:])
        ind = np.where(np.array(alpha_list) != 0.0)
        plt.figure()
        for j, mesh_density in enumerate(mesh_densities[:-1]):
            plt.plot(np.array(alpha_list)[ind], err[j,ind].flatten(), 'o', label=mesh_density.replace("_", " ").title(), markersize=10-3*j, mec='k', mfc='none')
        #for j, mesh_density in enumerate(mesh_densities):
        #    plt.plot(alpha_list, CL[i,j,:], 'o', label=mesh_density.replace("_", " ").title(), markersize=10-2*j, mec='k', mfc='none')
        #if case=="MH":
        #    plt.legend()
        plt.xlabel('$\\alpha\,[^\circ]$')
        plt.ylabel('Fractional Error in $C_y$')
        plt.yscale('log')
        plt.ylim(bottom=2e-4, top=1e-1)
        plot_loc = study_dir + 'plots/delta_wing_Cy_convergence{0}'.format(quad_label)
        plt.savefig(plot_loc + ".pdf")
        plt.savefig(plot_loc + ".svg")
        plt.close()
    
        # Plot most refined results against experimental data

        # CD plot
        plt.figure()
        plt.plot(alpha_list, CD[i,-1,:], 'ko', label="MachLine", markersize=3)
        plt.plot(CD_exp[:,0], CD_exp[:,1], 'ks', label='Experiment', markersize=3)
        plt.plot(CPanel_CD[:,0], CPanel_CD[:,1], 'kv', label='CPanel', markersize=3)
        if case=="MH":
            plt.legend()
        plt.xlabel('$\\alpha\,[^\circ]$')
        plt.ylabel('$C_D$')
        plt.ylim(bottom=0.0)
        plot_loc = study_dir + 'plots/delta_wing_CD_comparison_{0}{1}'.format(PRESSURE_FOR_FORCES, quad_label)
        plt.savefig(plot_loc + ".pdf")
        plt.savefig(plot_loc + ".svg")
        plt.close()

        # CL plot
        plt.figure()
        plt.plot(alpha_list, CL[i,-1,:], 'ko', label="MachLine", markersize=3)
        plt.plot(CL_exp[:,0], CL_exp[:,1], 'ks', label='Experiment', markersize=3)
        plt.plot(CPanel_CL[:,0], CPanel_CL[:,1], 'kv', label='CPanel', markersize=3)
        if case=="MH":
            plt.legend()
        plt.xlabel('$\\alpha\,[^\circ]$')
        plt.ylabel('$C_L$')
        plot_loc = study_dir + 'plots/delta_wing_CL_comparison_{0}{1}'.format(PRESSURE_FOR_FORCES, quad_label)
        plt.savefig(plot_loc + ".pdf")
        plt.savefig(plot_loc + ".svg")
        plt.close()

    # Analyze force convergence
    plt.figure()
    print()
    print("CL Convergence Rate")
    print("-------------------")
    for i, case in enumerate(cases):
        slopes = []
        alpha_to_plot = []
        for k, alpha in enumerate(alpha_list):
            if alpha==0.0:
                continue
            slopes.append(get_order_of_convergence(l_avg[1:], CL[i,1:,k], truth_from_results=True))
            alpha_to_plot.append(alpha)
        print("{0}: {1} +/- {2}".format(case, round(np.average(slopes), 4), round(np.std(slopes), 4)))
        plt.plot(alpha_to_plot, slopes, line_styles[i], label=case)
    plt.xlabel('$\\alpha$')
    plt.ylabel('Order of Convergence in $C_L$')
    plt.legend()
    plt.savefig(study_dir + "plots/CL_convergence_over_alpha.pdf")
    plt.savefig(study_dir + "plots/CL_convergence_over_alpha.svg")
    plt.close()

    plt.figure()
    print()
    print("CD Convergence Rate")
    print("-------------------")
    for i, case in enumerate(cases):
        slopes = []
        for k, alpha in enumerate(alpha_list):
            slopes.append(get_order_of_convergence(l_avg[1:], CD[i,1:,k], truth_from_results=True))
        print("{0}: {1} +/- {2}".format(case, round(np.average(slopes), 4), round(np.std(slopes), 4)))
        plt.plot(alpha_list, slopes, line_styles[i], label=case)
    plt.xlabel('$\\alpha$')
    plt.ylabel('Order of Convergence in $C_D$')
    plt.legend()
    plt.savefig(study_dir + "plots/CD_convergence_over_alpha.pdf")
    plt.savefig(study_dir + "plots/CD_convergence_over_alpha.svg")
    plt.close()

    plt.figure()
    print()
    print("Cx Convergence Rate")
    print("-------------------")
    for i, case in enumerate(cases):
        slopes = []
        for k, alpha in enumerate(alpha_list):
            slopes.append(get_order_of_convergence(l_avg[1:], Cx[i,1:,k], truth_from_results=True))
        print("{0}: {1} +/- {2}".format(case, round(np.average(slopes), 4), round(np.std(slopes), 4)))
        plt.plot(alpha_list, slopes, line_styles[i], label=case)
    plt.xlabel('$\\alpha$')
    plt.ylabel('Order of Convergence in $C_x$')
    plt.legend()
    plt.savefig(study_dir + "plots/Cx_convergence_over_alpha.pdf")
    plt.savefig(study_dir + "plots/Cx_convergence_over_alpha.svg")
    plt.close()

    plt.figure()
    print()
    print("Cy Convergence Rate")
    print("-------------------")
    for i, case in enumerate(cases):
        slopes = []
        alpha_to_plot = []
        for k, alpha in enumerate(alpha_list):
            if alpha==0.0:
                continue
            slopes.append(get_order_of_convergence(l_avg[1:], Cy[i,1:,k], truth_from_results=True))
            alpha_to_plot.append(alpha)
        print("{0}: {1} +/- {2}".format(case, round(np.average(slopes), 4), round(np.std(slopes), 4)))
        plt.plot(alpha_to_plot, slopes, line_styles[i], label=case)
    plt.xlabel('$\\alpha$')
    plt.ylabel('Order of Convergence in $C_z$')
    #plt.legend()
    plt.savefig(study_dir + "plots/Cy_convergence_over_alpha.pdf")
    plt.savefig(study_dir + "plots/Cy_convergence_over_alpha.svg")
    plt.close()


def plot_time_accuracy_comp(alpha_list):

    # Initialize storage
    Cx = np.zeros((len(cases), len(mesh_densities), len(alpha_list)))
    Cy = np.zeros((len(cases), len(mesh_densities), len(alpha_list)))
    t = np.zeros((len(cases), len(mesh_densities), len(alpha_list)))
    l_avg = np.zeros(len(mesh_densities))

    # Iterate over cases
    for i, (case, quad_label) in enumerate(zip(cases, quad_labels)):

        # Iterate over mesh densities
        for j, mesh_density in enumerate(mesh_densities):

            # Iterate over alpha
            for k, alpha in enumerate(alpha_list):

                # Get report
                report_loc = reports_dir + 'delta_wing_{0}_{1}{2}.json'.format(alpha, mesh_density, quad_label)
                with open(report_loc, 'r') as results_handle:
                    report = json.load(results_handle) 

                # Store force data
                force_coefs = report['total_forces']
                Cx[i,j,k] = force_coefs['Cx']
                Cy[i,j,k] = force_coefs['Cy']

                # Store average characteristic length
                l_avg[j] = report["mesh_info"]["average_characteristic_length"]

                # Store time
                t[i,j,k] = report["total_runtime"]

    # Calculate error
    err_Cx = np.abs((Cx[:,:-1,:] - Cx[:,-1,:])/Cx[:,-1,:])
    err_Cy = np.abs((Cy[:,:-1,:] - Cy[:,-1,:])/Cy[:,-1,:])

    for k, alpha in enumerate(alpha_list):

        # Plot Cx
        plt.figure()
        for i, (case, quad_label) in enumerate(zip(cases, quad_labels)):
            plt.plot(t[i,:-1,k], err_Cx[i,:,k], line_styles[i], label=case)
        plt.xscale('log')
        plt.yscale('log')
        if alpha == 0.0:
            plt.legend()
        plt.xlabel("Computation Time $[s]$")
        plt.ylabel("Fractional Error in $C_x$")
        plt.savefig(plot_dir + "Cx_accuracy_vs_time_aoa_{0}.pdf".format(int(alpha)))
        plt.savefig(plot_dir + "Cx_accuracy_vs_time_aoa_{0}.svg".format(int(alpha)))
        plt.close()

        # Plot Cy
        plt.figure()
        for i, (case, quad_label) in enumerate(zip(cases, quad_labels)):
            if alpha != 0.0:
                plt.plot(t[i,:-1,k], err_Cy[i,:,k], line_styles[i], label=case)
            else:
                plt.plot(t[i,:,k], abs(Cy[i,:,k]), line_styles[i], label=case)
        plt.xscale('log')
        plt.yscale('log')
        if alpha == 0.0:
            plt.legend()
        plt.xlabel("Computation Time $[s]$")
        if alpha != 0.0:
            plt.ylabel("Fractional Error in $C_z$")
        else:
            plt.ylabel("Absolute Error in $C_z$")
        plt.savefig(plot_dir + "Cy_accuracy_vs_time_aoa_{0}.pdf".format(int(alpha)))
        plt.savefig(plot_dir + "Cy_accuracy_vs_time_aoa_{0}.svg".format(int(alpha)))
        plt.close()


def run_pressure_distribution_comparison(mesh_density):

    # Parameters
    angles_of_attack = [0.0, 2.0, 4.1, 8.5, 10.75]

    # Run MachLine
    run_machline_cases(angles_of_attack, M, mesh_density)

    # Plot pressure rule method over a range of angles of attack at each semispan location
    pressure_rule = "isentropic" # isentropic, second-order, slender-body, or linear
    plot_pressure_slices(pressure_rule, angles_of_attack, mesh_density)


def run_force_comparison():

    # Parameters
    angles_of_attack = [-6., -5., -4., -3., -2., -1., 0., 1., 2., 3., 4., 5.]

    # Run MachLine
    #for mesh_density in mesh_densities:
    #    run_machline_cases(angles_of_attack, M, mesh_density)

    # Plot
    plot_force_convergence_over_alpha(angles_of_attack)

    # Plot time/accuracy comparison
    plot_time_accuracy_comp(angles_of_attack)


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
                "pressure_for_forces" : PRESSURE_FOR_FORCES
            },
            "output" : {
                "body_file" :          body_file,
                "control_point_file" : results_dir + "delta_wing_{0}_deg_{1}_control_points.vtk".format(alpha, mesh_density),
                "wake_file": results_dir + "delta_wing_{0}_deg_{1}_wake.vtk".format(alpha, mesh_density),
                "report_file": reports_dir + "delta_wing_{0}_{1}.json".format(alpha, mesh_density)
            }
        }

        # Dump
        input_file = study_dir + "input_files/delta_wing_input.json"
        write_input_file(input_dict, input_file)

        # Run
        run_quad(input_file, run=RERUN_MACHLINE)
        

if __name__=="__main__":

    # Run comparisons
    #for mesh_density in mesh_densities:
    #    run_pressure_distribution_comparison(mesh_density)
    run_force_comparison()