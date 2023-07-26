import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0,'studies/')
sys.path.insert(1,'dev/helper_scripts/')
sys.path.insert(2,'C:/Program Files/ParaView 5.11.0/bin/Lib/site-packages')

from case_running_functions import *
from off_body_points_creator import *
from paraview_functions import get_data_from_csv, get_data_column_from_array

RERUN_MACHLINE = True

study_dir = 'studies/supersonic_cone_velocities/'

result_dir = study_dir + 'results/'
report_dir = study_dir + 'reports/'
plot_dir = study_dir + 'plots/'
meshes = study_dir + 'meshes/'


def create_points(ray_angles, n_points, cone_angle):

    tip = [1/np.tan(np.deg2rad(cone_angle)),0.0, 0.0]

    for i, ray_angle in enumerate(ray_angles):
        # Create Points
        points = cone_ray_slice(tip, ray_angle, n_points)

        # Write points to csv
        points_file = study_dir + 'cone_{0}_deg_{1}_deg_ray_angle_points.csv'.format(cone_angle,ray_angle)
        write_points(points_file, points)
    


def run(grid, ray_angle, cone_angle):
    # Create input file
    case_name = "cone_{0}_deg_{1}".format(cone_angle,grid)

    result_file = result_dir + case_name + ".vtk"
    report_file = report_dir + case_name + ".json"
    points_file = study_dir + 'cone_{0}_deg_{1}_deg_ray_angle_points.csv'.format(cone_angle,ray_angle)
    points_output_file = result_dir + case_name + "_{0}_deg_ray_angle_points.csv".format(ray_angle)

    input_dict = {
        "flow" : {
            "freestream_velocity" : [-1.0,0.0,0.0],
            "freestream_mach_number" : 1.2
        },
        "geometry" : {
            "file" : meshes + case_name + "_improved.vtk",
            "spanwise_axis" : "+z",
            "mirror_about" : "xy"
        },
        "solver":{
        },
        "post_processing" : {
        },
        "output" : {
            "body_file" : result_file,
            "report_file" : report_file,
            "offbody_points" : {
                "points_file" : points_file,
                "output_file" : points_output_file
            }
        }
    }

    input_filename = study_dir + "input.json"
    write_input_file(input_dict,input_filename)
    reports = run_quad(input_filename, run=RERUN_MACHLINE)

    return reports


def plot_velocities(grids, ray_angles, cone_angle):
    
    line_colors = ['#BBBBBB', '#777777', '#000000']
    line_styles_local = ['-', '--', '-.']

    for k, (case, quad_label) in enumerate(zip(cases, quad_labels)):

        plt.figure()
        for i, (ray_angle, style) in enumerate(zip(ray_angles, line_styles_local)):
            for j, grid in enumerate(grids):
                filename = result_dir + "cone_{0}_deg_{1}_{2}_deg_ray_angle_points{3}.csv".format(cone_angle, grid, ray_angle, quad_label)

                headers, data = get_data_from_csv(filename, False)
                x = get_data_column_from_array(headers,data,"x")
                x = x/max(x)
                V_x = get_data_column_from_array(headers,data,"V_x")
                V_y = get_data_column_from_array(headers,data,"V_y")
                V_z = get_data_column_from_array(headers,data,"V_z")

                V_mag = np.zeros(len(V_x))
                for l in range(len(V_x)):
                    V_mag[l] = np.sqrt(float(V_x[l])**2 + float(V_y[l])**2 + float(V_z[l])**2)

                if j == 2:
                    plt.plot(x, V_mag, style, label='{0} deg'.format(ray_angle), color=line_colors[j], linewidth=1)
                else:
                    plt.plot(x, V_mag, style, color=line_colors[j], linewidth=1)

        plt.xlabel('$x/l$')
        plt.ylabel('Velocity Magnitude')
        plt.legend(title="Ray Angle", fontsize=6, title_fontsize=6)
        
        plt.savefig(plot_dir+"supersonic_cone_{0}_deg_velocities_{1}_mir.pdf".format(cone_angle, case))

        plt.close()


if __name__ == "__main__":

    # Make directories
    if (not os.path.exists(result_dir)):
        os.mkdir(result_dir)
    if (not os.path.exists(report_dir)):
        os.mkdir(report_dir)
    if (not os.path.exists(plot_dir)):
        os.mkdir(plot_dir)

    n_points = 50
    grids = ["coarse","medium","fine"]
    cone_angle = 5
    ray_angles = [10,15,20]

    create_points(ray_angles, n_points, cone_angle)

    for i, ray_angle in enumerate(ray_angles):
        for j, grid in enumerate(grids):
            run(grid, ray_angle, cone_angle)
        
    plot_velocities(grids, ray_angles, cone_angle)
    
