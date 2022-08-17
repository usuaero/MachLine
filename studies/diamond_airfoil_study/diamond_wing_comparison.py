import json
import subprocess as sp
import flow54 as fl
import numpy as np
import paraview.simple as pvs
import matplotlib.pyplot as plt


def run_comparison(M, alpha, grid, half_angle, run_machline=True):
    # Runs the comparison of the diamond wing to shock-expansion theory

    # Parameters
    R_G = 287.058
    gamma = 1.4
    T_inf = 300.0
    c_inf = np.sqrt(gamma*R_G*T_inf)
    rho = 1.225
    p_inf = 1.0e5

    # Storage locations
    case_name = "M_{0}_aoa_{1}_{2}_deg_{3}".format(M, alpha, int(half_angle), grid)
    plot_dir = "studies/diamond_airfoil_study/plots/"
    mesh_file = "studies/diamond_airfoil_study/meshes/diamond_{0}_deg_full_{1}.stl".format(int(half_angle), grid)
    results_file = "studies/diamond_airfoil_study/results/"+case_name+".vtk"
    report_file = "studies/diamond_airfoil_study/reports/"+case_name+".json"
    data_file = 'studies/diamond_airfoil_study/data/'+case_name+'.csv'

    if run_machline:

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
                "formulation": "morino",
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
                "body_file" : results_file,
                "report_file" : report_file
            }
        }

        # Dump
        input_file = "studies/diamond_airfoil_study/diamond_input.json"
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
    CL_se, CD_se = airfoil.get_inviscid_lift_and_drag_coefs()
    
    # Load MachLine report file
    with open(report_file, 'r') as report_handle:
        report = json.load(report_handle)

    # Extract lift and drag coefficients
    Cx = report["total_forces"]["Cx"]
    Cz = report["total_forces"]["Cz"]
    C_a = np.cos(np.radians(alpha))
    S_a = np.sin(np.radians(alpha))
    CL_ml = Cz*C_a - Cx*S_a
    CD_ml = Cx*C_a + Cz*S_a

    # Read into ParaView
    data_reader = pvs.LegacyVTKReader(registrationName=case_name, FileNames=results_file)

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
    pvs.SaveData(data_file, proxy=plot, PointDataArrays=data_to_process, FieldAssociation='Point Data', Precision=12)

    # Read in data
    data = np.genfromtxt(data_file, delimiter=',', skip_header=1)
    
    # Plot data from MachLine
    plt.figure()
    plt.plot(data[:,4], data[:,0], 'ks', markersize=3, label='Second-Order')
    plt.plot(data[:,4], data[:,1], 'kv', markersize=3, label='Isentropic')
    plt.plot(data[:,4], data[:,2], 'ko', markersize=3, label='Slender-Body')
    plt.plot(data[:,4], data[:,3], 'k^', markersize=3, label='Linear')

    # Plot data from shock-expansion theory
    x = np.linspace(0.0, 1.0, 100)
    Cp_upper = np.ones_like(x)
    Cp_upper[:50] *= Cp3
    Cp_upper[50:] *= Cp5
    Cp_lower = np.ones_like(x)
    Cp_lower[:50] *= Cp2
    Cp_lower[50:] *= Cp4
    plt.plot(x, Cp_upper, 'k-', label='Shock-Expansion', linewidth=0.5)
    plt.plot(x, Cp_lower, 'k-', linewidth=0.5)

    # Format
    plt.xlabel('$x$')
    plt.ylabel('$C_p$')
    plt.legend(title='Pressure Rule', fontsize=6, title_fontsize=6)
    plt.savefig(plot_dir+case_name+".pdf".format(M, alpha, int(half_angle)))
    plt.close()

    return CL_se, CD_se, CL_ml, CD_ml


if __name__=="__main__":

    grids = ["coarse", "medium", "fine", "ultra_fine"]
    N_verts = [402, 1090, 3770, 13266]
    Ms = [1.5, 2.0, 3.0, 5.0]
    alphas = np.linspace(0.0, 5.0, 6)
    half_angles = [5]

    CLs = np.zeros((len(grids), len(Ms), len(alphas), len(half_angles)))
    CDs = np.zeros((len(grids), len(Ms), len(alphas), len(half_angles)))

    for i, grid in enumerate(grids):
        for j, M in enumerate(Ms):
            for k, alpha in enumerate(alphas):
                for l, half_angle in enumerate(half_angles):

                    _,_,CLs[i,j,k,l], CDs[i,j,k,l] = run_comparison(M, alpha, grid, half_angle, run_machline=True)

    plot_dir = "studies/diamond_airfoil_study/plots/"
    for j, M in enumerate(Ms):
        for k, alpha in enumerate(alphas):
            for l,half_angle in enumerate(half_angles):

                case_name = "M_{0}_aoa_{1}_{2}_deg_convergence".format(M, alpha, int(half_angle))

                plt.figure()
                plt.plot(N_verts[:-1], abs((CLs[:-1,j,k,l]-CLs[-1,j,k,l])/CLs[-1,j,k,l]), 'k-')
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel('$N_{verts}$')
                plt.ylabel('Error in $C_L$')
                plt.savefig(plot_dir+case_name+"_CL.pdf")
                plt.close()

                plt.figure()
                plt.plot(N_verts[:-1], abs((CDs[:-1,j,k,l]-CDs[-1,j,k,l])/CDs[-1,j,k,l]), 'k-')
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel('$N_{verts}$')
                plt.ylabel('Error in $C_D$')
                plt.savefig(plot_dir+case_name+"_CD.pdf")
                plt.close()