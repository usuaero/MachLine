import json
import numpy as np
import subprocess as sp
import flow54 as fl
import numpy as np
import paraview.simple as pvs
import matplotlib.pyplot as plt


def run_sphere_comparison(grid, sample, run_machline=True):
    # Runs the comparison of the random sphere

    # Storage locations
    case_name = "random_sphere_{0}_sample_{1}".format(grid, sample)
    plot_dir = "studies/panel_regularity_sphere_study/plots/"
    body_file = "studies/panel_regularity_sphere_study/meshes/"+case_name+".vtk"
    report_file = "studies/panel_regularity_sphere_study/reports/"+case_name+".json"
    data_file = 'studies/panel_regularity_sphere_study/data/'+case_name+'.csv'

    if run_machline:

        # Declare MachLine input
        input_dict = {
            "flow" : {
                "freestream_velocity": [10.0, 0.0, 0.0]
            },
            "geometry" : {
                "file" : body_file,
                "spanwise_axis" : "+y",
                "wake_model" : {
                    "wake_present" : False
                },
                "reference" : {
                    "area" : np.pi
                }
            },
            "solver": {
                "formulation": "morino",
                "control_point_offset": 1.1e-6,
                "matrix_solver" : "BJAC",
                "relaxation" : 0.9,
                "tolerance" : 1.1e-11,
                "block_size" : 400
            },
            "post_processing" : {
            },
            "output" : {
                "body_file" : body_file,
                "report_file" : report_file
            }
        }

        # Dump
        input_file = "studies/panel_regularity_sphere_study/sphere_input.json"
        with open(input_file, 'w') as input_handle:
            json.dump(input_dict, input_handle, indent=4)

        # Run
        sp.run(["./machline.exe", input_file])
    
    # Load MachLine report file
    with open(report_file, 'r') as report_handle:
        report = json.load(report_handle)

    # Extract lift and drag coefficients
    Cx = report["total_forces"]["Cx"]
    Cy = report["total_forces"]["Cy"]
    Cz = report["total_forces"]["Cz"]

    # Read into ParaView
    sphere_vtk = pvs.LegacyVTKReader(registrationName=case_name, FileNames=body_file)
    pvs.SaveData(data_file, proxy=sphere_vtk, FieldAssociation="Cell Data")

    ## Filter cell data to point data
    #filter = pvs.CellDatatoPointData(registrationName='Filter', Input=sphere_vtk)
    #data_to_process = ['C_p_inc']
    #filter.CellDataArraytoprocess = data_to_process

    ## Extract and save data
    #plot = pvs.PlotOnIntersectionCurves(registrationName='Plot', Input=filter)
    #plot.SliceType = 'Plane'
    #plot.SliceType.Normal = [0.0, 1.0, 0.0]
    #view = pvs.CreateView('XYChartView')
    #display = pvs.Show(plot, view, 'XYChartRepresentation')
    #view.Update()
    #display.XArrayName = 'Points_X'
    #view.Update()
    #pvs.SaveData(data_file, proxy=plot, PointDataArrays=data_to_process, FieldAssociation='Point Data', Precision=12)

    # Read in data
    data = np.genfromtxt(data_file, delimiter=',', skip_header=1)
    
    # Plot data from MachLine
    plt.figure()
    theta = np.degrees(np.arccos(data[:,4]))
    theta_a = np.linspace(0.0, np.pi, 100)
    plt.plot(theta, data[:,7], 'k.', markersize=1, label='Random Mesh')
    plt.plot(np.degrees(theta_a), 1.0-2.25*np.sin(theta_a)**2, 'k-', label='Analytic')
    plt.xlabel('$\\theta [^\circ]$')
    plt.ylabel('$C_p$')
    plt.ylim(bottom=-2.5)
    plt.legend(fontsize=6, title_fontsize=6)
    plt.savefig(plot_dir+case_name+".pdf")
    plt.close()

    return Cx, Cy, Cz


if __name__=="__main__":

    Ns = [125, 250, 500, 1000, 2000]
    grids = ["ultra_coarse", "coarse", "medium", "fine", "ultra_fine"]
    samples = [x for x in range(10)]

    C_f = np.zeros((len(grids), len(samples), 3))

    # Plot each sample
    for i, grid in enumerate(grids):
        for j, sample in enumerate(samples):

                C_f[i,j,0], C_f[i,j,1], C_f[i,j,2] = run_sphere_comparison(grid, sample, run_machline=False)

    # Plot convergence
    avg = np.average(C_f, axis=1)
    std_dev = np.std(C_f, axis=1)
    plot_dir = "studies/panel_regularity_sphere_study/plots/"
    plt.figure()
    plt.plot(Ns, np.abs(avg[:,0]), 'ks-', label="$C_x$")
    plt.plot(Ns, np.abs(avg[:,1]), 'ko--', label="$C_y$")
    plt.plot(Ns, np.abs(avg[:,2]), 'kv-.', label="$C_z$")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("$N_{verts}$")
    plt.ylabel("Force Coefficient")
    plt.legend()
    plt.savefig("studies/panel_regularity_sphere_study/plots/convergence.pdf")