import os
import json

import numpy as np
import subprocess as sp
import paraview.simple as pvs
import matplotlib.pyplot as plt


def run_machline(input_dict, input_filename):
    # Runs MachLine using the provided input dictionary. Will write out to the specified filename

    # Write input
    with open(input_filename, 'w') as input_handle:
        json.dump(input_dict, input_handle)

    # Run MachLine
    sp.run(["./machline.exe", input_filename])

    # Get report
    report_file = input_dict["output"].get("report_file")
    if report_file is not None:
        with open(report_file) as report_handle:
            report = json.load(report_handle)
            return report


def extract_pressures(result_file):

    # Read into ParaView
    sphere_vtk = pvs.LegacyVTKReader(registrationName='sphere', FileNames=result_file)
    pvs.SaveData("temp.csv", proxy=sphere_vtk, FieldAssociation="Cell Data")

    # Read in data
    data = np.genfromtxt("temp.csv", delimiter=',', skip_header=1)
    #os.remove("temp.csv")

    # Get locations and pressures
    locs = data[:,4:7]
    locs[:,0] -= 1.0
    C_P = data[:,7]

    return locs, C_P


if __name__=="__main__":

    # Options
    density = "medium"
    phis = np.radians([30.0, 45.0, 60.0])
    thetas = np.radians([30.0, 45.0, 60.0])

    for i, phi in enumerate(phis):
        for j, theta in enumerate(thetas):

            # Initialize input
            result_file = "studies/convergence_study/results/sphere_{0}.vtk".format(density)
            report_file = "studies/convergence_study/reports/sphere_{0}.json".format(density)

            V_inf = [np.cos(phi)*np.cos(theta), np.sin(phi)*np.cos(theta), np.sin(theta)]

            # Lower-order results
            input_dict ={
                "flow" : {
                    "freestream_velocity" : V_inf
                },
                "geometry" : {
                    "file" : "studies/convergence_study/meshes/sphere_{0}.stl".format(density),
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

            # Run
            report = run_machline(input_dict, "studies/convergence_study/input.json")

            # Get pressures
            locs, C_P = extract_pressures(result_file)

            # Figure out thetas
            thetas = np.arccos(np.einsum('ij,j->i', locs, V_inf))
            C_P_anl = 1.0 - 2.25*np.sin(thetas)**2
            err = abs((C_P - C_P_anl)/C_P_anl)
            plt.figure()
            plt.plot(np.degrees(thetas), err, 'k.', markersize=1)
            plt.yscale('log')
            plt.xlabel('$\\theta [^\circ]$')
            plt.ylabel('$|(C_P - C_{P_{anl}})/C_{P_{anl}}|$')
            plt.show()

    ## plot
    #plt.figure()
    #for i in range(len(phis)):
    #    for j in range(len(thetas)):
    #        plt.plot(N, np.abs(Cx[i,j,:]), 'k-')
    #        plt.plot(N, np.abs(Cy[i,j,:]), 'k-')
    #        plt.plot(N, np.abs(Cz[i,j,:]), 'k-')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.xlabel('$N_{verts}$')
    #plt.ylabel('$C_F$')
    #plt.show()