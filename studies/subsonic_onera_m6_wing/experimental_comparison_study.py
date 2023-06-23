import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, './studies')
#sys.path.insert(1, 'C:/Program Files/ParaView 5.11.0/bin/Lib/site-packages')
from paraview_functions import get_data_from_csv, get_data_column_from_array, extract_plane_slice
from case_running_functions import write_input_file, run_quad, cases, quad_labels


RERUN_MACHLINE = False
study_dir = "studies/subsonic_onera_m6_wing/"


def run_wing_quad_comparison(AoA, mach_num, correction, grid):
    # Runs the Onera M6 wing at the angle of attack

    # Storage Locations
    if correction:
        case_name = "m6_onera_AoA_{0}_mach_{1}_{2}_pg".format(AoA, mach_num, grid)
    else:
        case_name = "m6_onera_AoA_{0}_mach_{1}_{2}_direct".format(AoA, mach_num, grid)
    mesh_file = study_dir + "meshes/M6_onera_{0}.stl".format(grid)
    #mesh_file = study_dir + "meshes/M6_onera.vtk"
    body_file = study_dir + "results/"+case_name+".vtk"
    control_point_file = body_file.replace(".vtk", "_control_points.vtk")
    report_file = study_dir + "reports/"+case_name+".json"

    # Calculate Freestream Velocity
    AoA_rd = ( AoA * np.pi ) / 180
    freestream_velocity = [np.cos(AoA_rd), 0.0, np.sin(AoA_rd)]

    # Pressure corrections
    correction_mach_number = 0.0
    prandtl_glauert = False

    if correction:
        correction_mach_number = mach_num
        prandtl_glauert = True
        mach_num = 0.0
        

    # Declare Machline Input
    input_dict = {
        "flow" : {
            "freestream_velocity" : freestream_velocity,
            "freestream_mach_number" : mach_num
        },
        "geometry" : {
            "file" : mesh_file,
            "mirror_about" : "xz",
            "spanwise_axis" : "+y",
            "reference" : {
                "area" : 0.7532
            }
        },
        "solver": {
            "write_A_and_b" : True
        },
        "post_processing" : {
            "pressure_rules": {
            },
            "subsonic_pressure_correction" : {
                "correction_mach_number" : correction_mach_number,
                "prandtl-glauert" : prandtl_glauert
            },
        },
        "output": {
            "body_file": body_file,
            "control_point_file": control_point_file,
            "report_file": report_file
        }
    }

    # Add pressure rule
    if correction:
        input_dict["post_processing"]["pressure_rules"]["incompressible"] = True
    else:
        input_dict["post_processing"]["pressure_rules"]["isentropic"] = True

    # Dump
    input_file = study_dir + "M6_input.json"
    write_input_file(input_dict,input_file)

    # Run
    reports = run_quad(input_file, run = RERUN_MACHLINE)

    # Get times
    t = []
    for report in reports:
        t.append(report["total_runtime"])

    return case_name, t


def get_min_and_max_SV_from_A_mat():
    # Reads in the A matrix and returns its largest and smallest singular values

    # Read in
    filename = "A_mat.txt"
    A = np.genfromtxt(filename)
    
    # Get singular values
    S = np.linalg.svd(A, compute_uv=False)

    return S[0], S[-1]


def run_cases(AoA, mach_num, grid):

    # Run With and without prandtl-glauert correction
    case_name_pg, _ = run_wing_quad_comparison(AoA, mach_num, True, grid)
    case_name_direct, t_direct = run_wing_quad_comparison(AoA, mach_num, False, grid)
    #S_max, S_min = get_min_and_max_SV_from_A_mat()
    #cond = np.sqrt(S_max/S_min)

    # Define Plot dir 
    plot_dir = study_dir + "plots/"

    # Pull Experimental Data
    column_headers, cell_data = get_data_from_csv(csv_file = study_dir + "experimental_data/M6_Data_mach_{0:.4f}".format(mach_num)+"_AoA_"+str(AoA)+".csv", remove_csv =False)
    
    # Define locations of y_slices
    y_loc = [0.239,0.526,0.778,0.957,1.077,1.136,1.184]
    semispan = [20,44,65,80,90,95,99]

    ## Loop Through Slices
    #for i, y in enumerate(y_loc):

    #    # Define Slice Values
    #    slice_normal = [0,1,0]
    #    slice_origin = [0,y,0]
    #    
    #    # Retrieve Experimental Data
    #    x_exp = get_data_column_from_array(column_headers,cell_data,"x/l_{0}".format(i+1))
    #    Cp_exp = get_data_column_from_array(column_headers,cell_data,"Cp_{0}".format(i+1))

    #    # Loop Through Formulations
    #    for j, (case, label) in enumerate(zip(cases, quad_labels)):
    #        
    #        # Pull Machline Data at slice
    #        data_dir_pg = study_dir + "results/" + case_name_pg + label + ".vtk"
    #        data_dir_direct = study_dir + "results/" + case_name_direct + label + ".vtk"

    #        # Switch between filtering and not
    #        for filter, which in zip([True, False], ['point', 'cell']):

    #            # Get slice
    #            headers_pg, slice_data_pg = extract_plane_slice(data_dir_pg, slice_normal, slice_origin, filter=filter, which_data=which)
    #            headers_direct, slice_data_direct = extract_plane_slice(data_dir_direct, slice_normal, slice_origin, filter=filter, which_data=which)

    #            # Get data we want to plot
    #            if filter:
    #                x_mach_pg = get_data_column_from_array(headers_pg,slice_data_pg,"Points:0")
    #                x_mach_direct = get_data_column_from_array(headers_direct,slice_data_direct,"Points:0")
    #            else:
    #                x_mach_pg = get_data_column_from_array(headers_pg,slice_data_pg,"centroid:0")
    #                x_mach_direct = get_data_column_from_array(headers_direct,slice_data_direct,"centroid:0")
    #            Cp_pg = get_data_column_from_array(headers_pg,slice_data_pg, "C_p_PG")
    #            Cp_ise_direct = get_data_column_from_array(headers_direct,slice_data_direct, "C_p_ise")

    #            # Modify machline x-axis
    #            x_mach_pg = [(x - min(x_mach_pg)) for x in x_mach_pg]
    #            x_mach_pg = [x/max(x_mach_pg) for x in x_mach_pg]
    #        
    #            x_mach_direct = [(x - min(x_mach_direct)) for x in x_mach_direct]
    #            x_mach_direct = [x/max(x_mach_direct) for x in x_mach_direct]

    #            # Differentiate between upper and lower surface
    #            half_pg = round(len(x_mach_pg)/2)
    #            half_direct = round(len(x_mach_direct)/2)

    #            # Plot Experimental data
    #            plt.figure()
    #            plt.plot(x_exp, Cp_exp, 'ko', label='Experiment', markersize=3)

    #            # Plot Machline Direct Data
    #            plt.plot(x_mach_direct[:half_direct],Cp_ise_direct[:half_direct], 'k-', linewidth=0.5, label='Direct')
    #            plt.plot(x_mach_direct[half_direct:],Cp_ise_direct[half_direct:], 'k-', linewidth=0.5)

    #            # Plot Machline Prandtl Glauert Correction Data
    #            plt.plot(x_mach_pg[:half_pg],Cp_pg[:half_pg], 'k--', linewidth=0.5, label='P-G')
    #            plt.plot(x_mach_pg[half_pg:],Cp_pg[half_pg:], 'k--', linewidth=0.5)

    #            # Define Case Name
    #            case_name = case_name_pg[:-3]

    #            # Plot
    #            plt.xlabel('$x/c$')
    #            plt.ylabel('$C_P$')
    #            plt.gca().invert_yaxis()
    #            plt.ylim(top=np.nanmin(Cp_exp).item()-0.5, bottom=np.nanmax(Cp_exp).item()+0.5)
    #            if j==1:
    #                plt.legend(fontsize=6, title_fontsize=6)

    #            # Save
    #            if not os.path.exists(study_dir + "plots/filtered/"):
    #                os.makedirs(study_dir + "plots/filtered/")
    #            if not os.path.exists(study_dir + "plots/unfiltered/"):
    #                os.makedirs(study_dir + "plots/unfiltered/")
    #            if filter:
    #                plt.savefig(plot_dir+"filtered/"+case_name+"_"+case+"_filtered_{0}.pdf".format(semispan[i]))
    #                plt.savefig(plot_dir+"filtered/"+case_name+"_"+case+"_filtered_{0}.svg".format(semispan[i]))
    #            else:
    #                plt.savefig(plot_dir+"unfiltered/"+case_name+"_"+case+"_{0}.pdf".format(semispan[i]))
    #                plt.savefig(plot_dir+"unfiltered/"+case_name+"_"+case+"_{0}.svg".format(semispan[i]))
    #            plt.close()

    ##return S_max, S_min, cond
    return t_direct


if __name__ == "__main__":

    # Define Case Mach Numbers
    mach_nums = [0.6971, 0.6977, 0.6990, 0.7001, 0.7003, 0.7009, 0.7019, 0.8399]
    AoA_nums = [6.09, 0.06, 3.06, 2.06, 1.08, 4.08, 5.06, 0.04]
    grids = ["coarse", "medium", "fine"]
    #grids = ["fine"]

    # Run Machline for all Mach Numbers and associated AoA
    print()
    print("Grid   M   alpha   S_max    S_min   cond")
    t = np.zeros((len(grids), len(mach_nums), 4))
    for j, grid in enumerate(grids):
        for i in range(len(mach_nums)):
            #S_max, S_min, cond = run_cases(AoA_nums[i], mach_nums[i], grid)
            #print(grid, mach_nums[i], AoA_nums[i], S_max, S_min, cond)
            t[j,i,:] = run_cases(AoA_nums[i], mach_nums[i], grid)

    # Print timing
    t_avg = np.average(t, axis=1)
    for i, density in enumerate(grids):
        for j, case in enumerate(cases):
            print(density, case, t_avg[i,j])
