import os
import json

import numpy as np
import subprocess as sp
import paraview.simple as pvs
import matplotlib.pyplot as plt

from studies.case_running_functions import run_machline, write_input_file
from studies.paraview_functions import extract_all_data, get_data_column_from_array


RERUN_MACHLINE = True


def extract_pressures(result_file):

    # Get data
    headers, data = extract_all_data(result_file, which_data='cell')

    # Get locations and pressures
    C_P = get_data_column_from_array(headers, data, 'C_p_inc')
    locs = np.zeros((len(C_P),3))
    locs[:,0] = get_data_column_from_array(headers, data, 'centroid:0') - 1.0
    locs[:,1] = get_data_column_from_array(headers, data, 'centroid:1')
    locs[:,2] = get_data_column_from_array(headers, data, 'centroid:2')

    return locs, C_P


if __name__=="__main__":

    # Options
    density = "medium"
    phis = np.radians([30.0, 45.0, 60.0])
    thetas = np.radians([30.0, 45.0, 60.0])

    for i, phi in enumerate(phis):
        for j, theta in enumerate(thetas):

            # Initialize input
            result_file = "studies/sphere/results/sphere_{0}.vtk".format(density)
            report_file = "studies/sphere/reports/sphere_{0}.json".format(density)

            V_inf = [np.cos(phi)*np.cos(theta), np.sin(phi)*np.cos(theta), np.sin(theta)]

            # Lower-order results
            input_dict ={
                "flow" : {
                    "freestream_velocity" : V_inf
                },
                "geometry" : {
                    "file" : "studies/sphere/meshes/sphere_{0}.stl".format(density),
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
            input_filename = "studies/sphere/input.json"
            write_input_file(input_dict, input_filename)
            report = run_machline(input_filename, run=RERUN_MACHLINE, delete_input=False)

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