import numpy as np
import matplotlib.pyplot as plt

from studies.case_running_functions import run_quad, get_order_of_convergence, cases, line_styles, write_input_file


RERUN_MACHLINE = False
study_dir = "studies/sphere/"


def run_cases_for_orientation_and_mesh_density(psi, theta, density):

    # Initialize input
    result_file = study_dir + "results/sphere_{0}.vtk".format(density)
    report_file = study_dir + "reports/sphere_{0}.json".format(density)

    # Lower-order results
    input_dict ={
        "flow" : {
            "freestream_velocity" : [np.cos(psi)*np.cos(theta), np.sin(psi)*np.cos(theta), np.sin(theta)]
        },
        "geometry" : {
            "file" : study_dir + "meshes/sphere_{0}.stl".format(density),
            "spanwise_axis" : "+y"
        },
        "solver" : {
        },
        "post_processing" : {
        },
        "output" : {
            "body_file" : result_file,
            "report_file" : report_file
        }
    }

    # Write to file
    write_input_file(input_dict, input_filename)

    # Run quad
    reports = run_quad(input_filename, run=RERUN_MACHLINE)

    N_sys = reports[0]["solver_results"]["system_dimension"]
    l_avg = reports[0]["mesh_info"]["average_characteristic_length"]

    C_F = np.zeros((4,3))
    for i, report in enumerate(reports):
        C_F[i,0] = report["total_forces"]["Cx"]
        C_F[i,1] = report["total_forces"]["Cy"]
        C_F[i,2] = report["total_forces"]["Cz"]

    return N_sys, l_avg, C_F


if __name__=="__main__":

    # Options
    input_filename = "studies/sphere/input.json"
    densities = ["ultra_coarse", "very_coarse", "coarse", "medium"]
    psis = np.radians([30.0, 45.0, 60.0])
    thetas = np.radians([30.0, 45.0, 60.0])
    N = []
    l_avg = []
    Cx = np.zeros((len(psis), len(thetas), len(densities), 4))
    Cy = np.zeros((len(psis), len(thetas), len(densities), 4))
    Cz = np.zeros((len(psis), len(thetas), len(densities), 4))

    for i, psi in enumerate(psis):
        for j, theta in enumerate(thetas):
            for k, density in enumerate(densities):

                # Run
                N_sys_i, l_avg_i, C_F = run_cases_for_orientation_and_mesh_density(psi, theta, density)

                # Store
                if i==0 and j==0:
                    N.append(N_sys_i)
                    l_avg.append(l_avg_i)
                Cx[i,j,k,:] = C_F[:,0]
                Cy[i,j,k,:] = C_F[:,1]
                Cz[i,j,k,:] = C_F[:,2]

    # Calculate norms of force vectors
    C_F_norms = np.sqrt(Cx**2 + Cy**2 + Cz**2)

    # Calculate convergence rates
    orders = [[], [], [], []]
    order_Cx = [[], [], [], []]
    order_Cy = [[], [], [], []]
    order_Cz = [[], [], [], []]
    order_norm = [[], [], [], []]
    for i in range(len(psis)):
        for j in range(len(thetas)):
            for k in range(4):
                orders[k].append(get_order_of_convergence(l_avg, Cx[i,j,:,k], truth_from_results=False))
                orders[k].append(get_order_of_convergence(l_avg, Cy[i,j,:,k], truth_from_results=False))
                orders[k].append(get_order_of_convergence(l_avg, Cz[i,j,:,k], truth_from_results=False))
                order_Cx[k].append(get_order_of_convergence(l_avg, Cx[i,j,:,k], truth_from_results=False))
                order_Cy[k].append(get_order_of_convergence(l_avg, Cy[i,j,:,k], truth_from_results=False))
                order_Cz[k].append(get_order_of_convergence(l_avg, Cz[i,j,:,k], truth_from_results=False))
                order_norm[k].append(get_order_of_convergence(l_avg, C_F_norms[i,j,:,k], truth_from_results=False))
            
    # Create arrays
    orders = np.array(orders)
    order_Cx = np.array(order_Cx)
    order_Cy = np.array(order_Cy)
    order_Cz = np.array(order_Cz)
    order_norm = np.array(order_norm)

    # Report average orders of convergence
    avg_orders = np.average(orders, axis=1)
    avg_order_Cx = np.average(order_Cx, axis=1)
    avg_order_Cy = np.average(order_Cy, axis=1)
    avg_order_Cz = np.average(order_Cz, axis=1)
    avg_order_norm = np.average(order_norm, axis=1)
    print()
    print("Average Orders of Convergence")
    print("-----------------------------")
    print("Between Forces")
    for k, case in enumerate(cases):
        print(case, ": ", avg_orders[k], " +/- ", np.std(orders[k]))
    print("Cx")
    for k, case in enumerate(cases):
        print(case, ": ", avg_order_Cx[k], " +/- ", np.std(order_Cx[k]))
    print("Cy")
    for k, case in enumerate(cases):
        print(case, ": ", avg_order_Cy[k], " +/- ", np.std(order_Cy[k]))
    print("Cz")
    for k, case in enumerate(cases):
        print(case, ": ", avg_order_Cz[k], " +/- ", np.std(order_Cz[k]))
    print("Norm")
    for k, case in enumerate(cases):
        print(case, ": ", avg_order_norm[k], " +/- ", np.std(order_norm[k]))

    # Plot
    plt.figure()
    for i in range(len(psis)):
        for j in range(len(thetas)):
            for k, (line_style, case) in enumerate(zip(line_styles, cases)):
                if i==0 and j==0:
                    plt.plot(l_avg, C_F_norms[i,j,:,k], line_style, label=case)
                    #plt.plot(l_avg, np.abs(Cx[i,j,:,k]), line_style, label=case)
                else:
                    plt.plot(l_avg, C_F_norms[i,j,:,k], line_style)
                    #plt.plot(l_avg, np.abs(Cx[i,j,:,k]), line_style)
                #plt.plot(l_avg, np.abs(Cy[i,j,:,k]), line_style)
                #plt.plot(l_avg, np.abs(Cz[i,j,:,k]), line_style)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$l_{avg}$')
    plt.ylabel('$||C_F||$')
    plt.legend()
    plt.savefig(study_dir + "plots/convergence.pdf")
    plt.savefig(study_dir + "plots/convergence.svg")
    plt.close()