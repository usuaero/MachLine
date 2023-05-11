import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, './studies')
sys.path.insert(1, 'C:/Program Files/ParaView 5.11.0/bin/Lib/site-packages')
from paraview_functions import get_data_from_csv, get_data_column_from_array, extract_plane_slice
from case_running_functions import write_input_file, run_quad, cases

RERUN_MACHLINE = True



def run_wing_quad_comparison(AoA, mach_num, correction=False):
    # Runs the Onera M6 wing at the angle of attack

    # Storage Locations
    if correction: case_name = "m6_onera_AoA_{0}_mach_{1}_pg".format(AoA, mach_num)
    else : case_name = "m6_onera_AoA_{0}_mach_{1}_direct".format(AoA, mach_num)
    mesh_file = "studies/subsonic_onera_m6_wing/meshes/m6_onera.vtk"
    body_file = "studies/subsonic_onera_m6_wing/results/"+case_name+".vtk"
    report_file = "studies/subsonic_onera_m6_wing/results/"+case_name+".json"


    # Calculate Freestream Velocity

    AoA_rd = ( AoA * np.pi ) / 180
    freestream_velocity = [np.cos(AoA_rd), 0.0, np.sin(AoA_rd)]

    # Pressure corrections
    correction_mach_number = 0.0
    prandtl_glauert = False

    if correction:
        correction_mach_number = mach_num
        prandtl_glauert = True
        mach_num = 0
        

    # Declare Machline Input
    input_dict = {
        "flow" : {
            "freestream_velocity" : freestream_velocity,
            "freestream_mach_number" : mach_num
        },
        "geometry" : {
            "file" : mesh_file,
            "mirror_about" : "xz",
            "spanwise_axis" : "+y"
        },
        "solver": {
            "formulation" : "morino"
        },
        "post_processing" : {
            "pressure_for_forces" : "isentropic",
            "pressure_rules": {
                "incompressible": False,
                "isentropic": True,
                "second-order": True
            },
            "subsonic_pressure_correction" : {
                "correction_mach_number" : correction_mach_number,
                "prandtl-glauert" : prandtl_glauert
            },
        },
        "output": {
            "body_file": body_file,
            "report_file": report_file
        }

    }

    # Dump
    input_file = "studies/subsonic_onera_m6_wing/Input_files/m6_input.json"
    write_input_file(input_dict,input_file)

    # Run
    reports = run_quad(input_file, run = RERUN_MACHLINE)





if __name__ == "__main__":

    # Define Case Mach Numbers
    mach_nums = [0.6971, 0.6977, 0.6990, 0.7001, 0.7003, 0.7009, 0.7019, 0.8399]
    AoA_nums = [6.09, 0.06, 3.06, 2.06, 1.08, 4.08, 5.06, 0.04]

    # Run Machline for all Mach Numbers and associated AoA
    for i in range(len(mach_nums)):
        print ("Run ", AoA_nums[i], mach_nums[i])
    
    #run_wing_quad_comparison(mach_num=0.6971, AoA=6.09)

    # Pull Experimental Data
    column_headers, cell_data = get_data_from_csv(csv_file = "studies/subsonic_onera_m6_wing/Experimental_Data/M6_Data_mach_0.6971_AoA_6.09.csv", remove_csv =False)
    x_exp = get_data_column_from_array(column_headers,cell_data,"x/l_1")
    Cp_exp = get_data_column_from_array(column_headers,cell_data,"Cp_1")

    # Define Slice Values
    slice_normal = [0,1,0]
    slice_origin = [0,0.239,0]
    data_dir = "studies/subsonic_onera_m6_wing/results/m6_onera_AoA_6.09_mach_0.6971_direct_QUAD_higher-order_morino.vtk"

    # Pull Machline Data
    column_headers, slice_data = extract_plane_slice(data_dir,slice_normal,slice_origin,filter=True)
    x_mach = get_data_column_from_array(column_headers,slice_data,"centroid:0")
    Cp_ise_direct = get_data_column_from_array(column_headers,slice_data, "C_p_ise")

    # Modify machline x-axis
    x_mach = [(x - min(x_mach)) for x in x_mach]
    x_mach = [x/x_mach[-1 : ] for x in x_mach]

    # Differentiate between upper and lower surfaces of Machline results
    half = round(len(x_mach)/2)
    print(half)


    # Plot Experimental Data
    plot_dir = "studies/subsonic_onera_m6_wing/plots/"

    plt.figure()
    plt.plot(x_exp,Cp_exp, 'ks', label='Experimental')
    plt.plot(x_mach[0:half],Cp_ise_direct, 'k-', label='Direct')
    plt.plot(x_mach[half:],Cp_ise_direct, 'k-', label='Direct')
    plt.xlabel('$x/l$')
    plt.ylabel('$C_p$')
    plt.gca().invert_yaxis()
    plt.legend(fontsize=6, title_fontsize=6)
    plt.show()



