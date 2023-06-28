import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0,'studies/')
sys.path.insert(1,'dev/helper_scripts/')
sys.path.insert(2, 'C:/Program Files/ParaView 5.11.0/bin/Lib/site-packages')

from case_running_functions import *
from off_body_points_creator import *
from paraview_functions import get_data_from_csv, get_data_column_from_array, extract_plane_slice

RERUN_MACHLINE = True

study_dir = 'studies/sphere/'

results = study_dir + 'results/'
report_dir = study_dir + 'reports/'
plots = study_dir + 'plots/'
meshes = study_dir + 'meshes/'
    
points_file = study_dir + 'sphere_surface_slice_points.csv'

def create_points(radius, offset, N_points):
   # Create points
    points, thetas = sphere_surface_slice([1.0,0,0],radius,offset,NPoints)

    # Write points to csv
    write_points(points_file, points)

    return points, thetas

def run(grid,offset):
    # Create input file
    result_file = results + "sphere_" +grid+".vtk"
    report_file = report_dir + "sphere_"+grid+".json"
    points_output_file = results + "sphere_offbody_points_{0:4.0e}_{1}.csv".format(offset,grid)
    
    input_dict = {
        "flow" : {
            "freestream_velocity" : [-1.0,0.0,0.0],
            "freestream_mach_number" : 0.0
        },
        "geometry" : {
            "file" : study_dir + "meshes/sphere_"+grid+".stl",
            "spanwise_axis" : "+y"
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

    # Run
    input_filename = study_dir + "input.json"
    write_input_file(input_dict,input_filename)
    reports = run_quad(input_filename, run=RERUN_MACHLINE)

    return reports

def get_data(thetas, radius, offset, grid):
    V_exact = []
    # Calculate predicted velocity on points
    for i, theta in enumerate(thetas):
        mu = -2*np.pi*radius**3
        V_r = mu*np.cos(theta)/(2*np.pi*(radius+offset)**3)
        V_theta = mu*np.sin(theta)/(4*np.pi*(radius+offset)**3)
        V_x = 1 + V_r*np.cos(theta) - V_theta*np.sin(theta)
        V_y = 2*V_r*np.sin(theta) - V_theta*np.cos(theta)
        V_mag = np.sqrt(V_x**2 + V_y**2)
        V_exact.append(V_mag)

    V_exact = np.array(V_exact)
    #headers, data = extract_plane_slice(results+"sphere_"+grid+"_QUAD_higher-order_morino.vtk",[0,0,1],[1,0,0],filter=True,which_data="point")

    #surface_x = get_data_column_from_array(headers,data, "Points:0")
    #surface_y = get_data_column_from_array(headers,data, "Points:1")
    #surface_thetas = 2*np.arctan2(surface_x,surface_y)
    #V_x = get_data_column_from_array(headers,data, "v:0")
    #V_y = get_data_column_from_array(headers,data, "v:1")
    #V_z = get_data_column_from_array(headers,data, "v:2")
    #V_surface = np.sqrt(V_x**2 + V_y**2 + V_z**2)

    #V_surface = np.array(V_surface)


    # Retrieve data from MachLine
    V_MachLine =[[],[],[],[]]


    for i, (case,quad_label) in enumerate(zip(cases, quad_labels)):
        filename = results + "sphere_offbody_points_{0:4.0e}_".format(offset) + grid + quad_label + ".csv"

        headers, data = get_data_from_csv(filename, False)
        V_x = get_data_column_from_array(headers,data,"V_x")
        V_y = get_data_column_from_array(headers,data,"V_y")
        V_z = get_data_column_from_array(headers,data,"V_z")

        V_mag = np.zeros(len(V_x))
        for j in range(len(V_x)):
            V_mag[j] = np.sqrt(float(V_x[j])**2 + float(V_y[j])**2 + float(V_z[j])**2)

        V_MachLine[i] = V_mag

    return V_exact, V_MachLine

def get_avg_error(thetas, exact, approx, type="Fractional"):

    error = 0.0

    for i, theta in enumerate(thetas):
        temp_error = abs(exact[i] - approx[i])
        if type == "Fractional":
            temp_error = temp_error/exact[i]
        error += temp_error
    
    error = error / len(thetas)

    return error


def plot_velocities(thetas, V_exact, V_MachLine, offset,grid):
    # Plot Velocities
    # Plot Analytical Data
    plt.figure()
    plt.plot(thetas,V_exact, label='Exact')

    # Plot MachLine Data
    for i, (case,style) in enumerate(zip(cases,line_styles)):
        plt.plot(thetas,V_MachLine[i],style,label=case,linewidth=0.5)

    # Other plotting stuff
    plt.xlabel('$\\theta$')
    plt.ylabel('$V$')
    plt.legend(fontsize=6,title_fontsize=6)

    # Save
    plt.savefig(plots+"incompressible_sphere_velocities_offset_{0:4.0e}_{1}.pdf".format(offset,grid))
    plt.savefig(plots+"incompressible_sphere_velocities_offset_{0:4.0e}_{1}.svg".format(offset,grid))

    plt.close()

def plot_error(thetas, V_exact, V_MachLine, offset,grid,type="Fractional"):

    # Plot Error
    plt.figure()
    for i, (case,style) in enumerate(zip(cases,line_styles)):
        if type == "Fractional":
            plt.plot(thetas,abs(V_exact-V_MachLine[i])/V_exact,style,label=case,linewidth=0.5)
        elif type == "Absolute":
            plt.plot(thetas,abs(V_exact-V_MachLine[i]),style,label=case,linewidth=0.5)

    # Other Plotting
    plt.xlabel('$\\theta$')
    plt.ylabel(type + ' Error')
    plt.legend(fontsize=6,title_fontsize=6)
    plt.yscale('log')

    # Save
    plt.savefig(plots+"incompressible_sphere_"+type+"_error_offset_{0:4.0e}_{1}.pdf".format(offset,grid))
    plt.savefig(plots+"incompressible_sphere_"+type+"_error_offset_{0:4.0e}_{1}.svg".format(offset,grid))

    plt.close()

def plot_error_vs_offset(thetas,radius,offsets,grid):

    plt.figure()
    for j, (case,style) in enumerate(zip(cases,line_styles)):
        error = []
        for i, offset in enumerate(offsets):
            V_exact,V_Machline = get_data(thetas,radius,offset,grid)
            error.append(get_avg_error(thetas,V_exact,V_Machline[j]))
        
        plt.plot(offsets,error,style,label=case,linewidth=0.5)

    plt.xlabel('Offset')
    plt.ylabel('Avg. Fractional Error')
    plt.legend(fontsize=6, title_fontsize=6)
    plt.yscale('log')
    plt.xscale('log')

    # Save
    plt.savefig(plots+"incompressible_sphere_error_vs_offset_{0}.pdf".format(grid))
    plt.savefig(plots+"incompressible_sphere_error_vs_offset_{0}.svg".format(grid))

    plt.close()

def grid_convergence(offset,radius,grids):

    # Run
    lengths = []
    for i,grid in enumerate(grids):
        reports = run(grid,offset)
        lengths.append(reports[i]["mesh_info"]["average_characteristic_length"])

    # Plot
    plt.figure()
    for j,(case,style) in enumerate(zip(cases,line_styles)):
        error = []
        for i,grid in enumerate(grids):
            
            V_exact,V_machline = get_data(thetas,radius,offset,grid)
            error.append(get_avg_error(thetas,V_exact,V_machline[j]))

        plt.plot(lengths,error,style,label=case,linewidth=0.5)

    plt.xlabel('Avg. Characteristic Length')
    plt.ylabel('Avg. Fractional Error')
    plt.legend(fontsize=6, title_fontsize=6)
    plt.yscale('log')
    plt.xscale('log')

    # Save
    plt.savefig(plots+"incompressible_sphere_grid_convergence_{0:4.0e}.pdf".format(offset))
    plt.savefig(plots+"incompressible_sphere_grid_convergence_{0:4.0e}.svg".format(offset))
    





    

if __name__ == "__main__":

    # Make directories
    if (not os.path.exists(results)):
        os.mkdir(results)
    if (not os.path.exists(report_dir)):
        os.mkdir(report_dir)
    if (not os.path.exists(plots)):
        os.mkdir(plots)

    radius = 1
    offset = 1e-2
    NPoints = 1000
    grids = ["ultra_coarse","very_coarse","coarse","medium"]

    points, thetas = create_points(radius, offset, NPoints)
    #for i, grid in enumerate(grids):
        #run(grid,offset)
        #V_exact,V_MachLine = get_data(thetas, radius, offset,grid)
        #plot_velocities(thetas, V_exact,V_MachLine,offset,grid)
        #plot_error(thetas,V_exact, V_MachLine,offset,grid)
        #plot_error(thetas,V_exact,V_MachLine,offset,grid,type="Absolute")

    offsets = [1e-3,1e-2,1e-1,1e0]
    for i,grid in enumerate(grids):
        for j, offset in enumerate(offsets):
            points, thetas = create_points(radius,offset,NPoints)
            run(grid,offset)
        plot_error_vs_offset(thetas,radius,offsets,grid)
 
    #grid_convergence(offset,radius,grids)
