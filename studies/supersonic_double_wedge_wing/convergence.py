import json
import numpy as np
import matplotlib.pyplot as plt
from studies.case_running_functions import run_quad, get_order_of_convergence, write_input_file


RERUN_MACHLINE = False


def run_quad_for_mach_aoa_and_mesh(M, alpha, density):
    """Runs a case quad for the given Mach number, angle of attack, and mesh density."""

    # Parameters
    half_angle = 5
    R_G = 287.058
    gamma = 1.4
    T_inf = 300.0
    c_inf = np.sqrt(gamma*R_G*T_inf)
    rho = 1.225
    p_inf = 1.0e5

    # Storage locations
    case_name = "M_{0}_aoa_{1}_{2}_deg_{3}".format(M, alpha, int(half_angle), grid)
    plot_dir = "studies/supersonic_double_wedge_wing/plots/"
    mesh_file = "studies/supersonic_double_wedge_wing/meshes/diamond_{0}_deg_full_{1}.stl".format(int(half_angle), grid)
    results_file = "studies/supersonic_double_wedge_wing/results/"+case_name+".vtk"
    report_file = "studies/supersonic_double_wedge_wing/reports/"+case_name+".json"
    data_file = 'studies/supersonic_double_wedge_wing/data/'+case_name+'.csv'

    # Write out input file

    # Declare MachLine input
    input_dict = {
        "flow": {
            "freestream_velocity": [M*c_inf*np.cos(np.radians(alpha)), 0.0, M*c_inf*np.sin(np.radians(alpha))],
            "gamma" : gamma,
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "wake_model": {
                "append_wake" : False,
            },
            "reference": {
                "area": 4.0
            }
        },
        "solver": {
            "formulation": "morino"
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
            "body_file" : results_file,
            "report_file" : report_file
        }
    }

    # Dump
    input_file = "studies/supersonic_double_wedge_wing/input.json"
    write_input_file(input_dict, input_file)

    # Run quad
    reports = run_quad(input_file, run=RERUN_MACHLINE)

    # Pull out forces
    C_F = np.zeros((4,3))
    for i, report in enumerate(reports):
        C_F[i,0] = report["total_forces"]["Cx"]
        C_F[i,1] = report["total_forces"]["Cy"]
        C_F[i,2] = report["total_forces"]["Cz"]

    # Get system dimension and average characteristic length
    N_sys = reports[0]["solver_results"]["system_dimension"]
    l_avg = reports[0]["mesh_info"]["average_characteristic_length"]

    return N_sys, l_avg, C_F


if __name__=="__main__":

    # Study parameters
    grids = ["coarse", "medium", "fine", "ultra_fine"]
    Ms = [1.5, 2.0, 3.0, 5.0]
    alphas = np.linspace(0.0, 5.0, 6)

    # Loop through parameters
    N_sys = np.zeros(len(grids))
    l_avg = np.zeros(len(grids))
    C_F = np.zeros((len(grids), len(Ms), len(alphas), 4, 3))
    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):
            for k, alpha in enumerate(alphas):

                N_sys[i], l_avg[i], C_F[i,j,k] = run_quad_for_mach_aoa_and_mesh(M, alpha, grid)

    # Set up plotting
    plot_dir = "studies/supersonic_double_wedge_wing/plots/"

    #Analyze convergence of Cx
    slopes = []
    for l, case in enumerate(['ML', 'MH', 'SL', 'SH']):
        slopes.append([])
        plt.figure()
        for j, M in enumerate(Ms):
            for k, alpha in enumerate(alphas):
                #err = abs((CLs[:-1,j,k,l]-CLs[-1,j,k,l])/CLs[-1,j,k,l])
                #plt.plot(N_verts[:-1], err, 'k-', linewidth=1)
                #coefs = np.polyfit(np.log(N_verts[:-1]), np.log(err), deg=1)
                slopes[l].append(get_order_of_convergence(l_avg, C_F[:,j,k,l,0], truth_from_results=True))

        print("Average Cx convergence rate case {0}: ".format(case), np.average(slopes[l]))

    #Analyze convergence of Cz
    slopes = []
    for l, case in enumerate(['ML', 'MH', 'SL', 'SH']):
        slopes.append([])
        plt.figure()
        for j, M in enumerate(Ms):
            for k, alpha in enumerate(alphas):
                #err = abs((CLs[:-1,j,k,l]-CLs[-1,j,k,l])/CLs[-1,j,k,l])
                #plt.plot(N_verts[:-1], err, 'k-', linewidth=1)
                #coefs = np.polyfit(np.log(N_verts[:-1]), np.log(err), deg=1)
                slopes[l].append(get_order_of_convergence(l_avg, C_F[:,j,k,l,2], truth_from_results=True))

        print("Average Cz convergence rate case {0}: ".format(case), np.average(slopes[l]))


        #plt.xscale('log')
        #plt.yscale('log')
        #plt.xlabel('$N_{verts}$')
        #plt.ylabel('Percent Error in $C_L$')
        #plt.savefig(plot_dir+"collected_CL_convergence_{0}.pdf".format(case))
        #plt.savefig(plot_dir+"collected_CL_convergence_{0}.svg".format(case))
        #plt.close()

        ## Plot combined CD convergences
        #slopes = []
        #plt.figure()
        #for j, M in enumerate(Ms):
        #    for k, alpha in enumerate(alphas):
        #        for l,half_angle in enumerate(half_angles):
        #            err = abs((CLs[:-1,j,k,l]-CLs[-1,j,k,l])/CLs[-1,j,k,l])
        #            plt.plot(N_verts[:-1], abs((CDs[:-1,j,k,l]-CDs[-1,j,k,l])/CDs[-1,j,k,l]), 'k-', linewidth=1)
        #            coefs = np.polyfit(np.log(N_verts[:-1]), np.log(err), deg=1)
        #            slopes[l].append(coefs[0])

        #print("Average CD convergence rate: ", np.average(slopes))
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.xlabel('$N_{verts}$')
        #plt.ylabel('Percent Error in $C_D$')
        #plt.savefig(plot_dir+"collected_CD_convergence_{0}.pdf".format(case))
        #plt.savefig(plot_dir+"collected_CD_convergence_{0}.svg".format(case))
        #plt.close()