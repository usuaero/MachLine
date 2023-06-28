import sys
import time
import numpy as np
sys.path.insert(0, '../panairwrapper')
import panairwrapper.panairwrapper as pnr
import matplotlib.pyplot as plt

from dev.helper_scripts.geometry_creator import _export_vtk
from studies.case_running_functions import run_quad, write_input_file, cases, line_styles, get_order_of_convergence
from studies.paraview_functions import extract_plane_slice, get_data_column_from_array, get_data_from_csv

# Geometry parameters (for Love's delta wing #11)
t_max = 0.08
c_t_max = 0.18
c_root = 0.23
b_semi = 0.463/2.0
A_ref = c_root*b_semi

# Run parameters
study_dir = "studies/pan_air_comparison/"
plot_dir = study_dir + "plots/"
pan_air_data_file = study_dir + "panair_files/panair_data.csv"
pan_air_output_file = study_dir + "panair_files/panair.out"


def create_networks(N_chordwise_fore, N_chordwise_aft, N_spanwise):
    # Creates the networks of vertices for PAN AIR to use

    networks = []
    network_names = []

    # Initialize fore network
    fore_network = np.zeros((N_chordwise_fore, N_spanwise, 3))

    # Get points at root
    x_root = np.linspace(0.0, c_t_max*c_root, N_chordwise_fore)
    y_root = np.zeros(N_chordwise_fore)
    z_root = np.linspace(0.0, 0.5*t_max*c_root, N_chordwise_fore)

    # Get points at tip
    x_tip = c_root
    y_tip = b_semi
    z_tip = 0.0

    # Interpolate
    for i in range(N_chordwise_fore):
        fore_network[i,:,0] = np.linspace(x_root[i], x_tip, N_spanwise)
        fore_network[i,:,1] = np.linspace(y_root[i], y_tip, N_spanwise)
        fore_network[i,:,2] = np.linspace(z_root[i], z_tip, N_spanwise)

    # Add
    networks.append(fore_network)
    network_names.append("ufr")

    # Add bottom fore network
    bottom_fore_network = np.copy(fore_network)
    bottom_fore_network = bottom_fore_network[::-1,:,:]
    bottom_fore_network[:,:,2] = -bottom_fore_network[:,:,2]
    networks.append(bottom_fore_network)
    network_names.append("bfr")

    # Initialize aft network
    aft_network = np.zeros((N_chordwise_aft, N_spanwise, 3))

    # Get points at root
    x_root = np.linspace(c_t_max*c_root, c_root, N_chordwise_aft)
    y_root = np.zeros(N_chordwise_aft)
    z_root = np.linspace(0.5*t_max*c_root, 0.0, N_chordwise_aft)

    # Interpolate
    for i in range(N_chordwise_aft):
        aft_network[i,:,0] = np.linspace(x_root[i], x_tip, N_spanwise)
        aft_network[i,:,1] = np.linspace(y_root[i], y_tip, N_spanwise)
        aft_network[i,:,2] = np.linspace(z_root[i], z_tip, N_spanwise)

    # Add
    networks.append(aft_network)
    network_names.append("uaf")

    # Add bottom aft network
    bottom_aft_network = np.copy(aft_network)
    bottom_aft_network = bottom_aft_network[::-1,:,:]
    bottom_aft_network[:,:,2] = -bottom_aft_network[:,:,2]
    networks.append(bottom_aft_network)
    network_names.append("baf")

    # Add wake network
    wake_network = np.zeros((2, N_spanwise, 3))
    wake_network[0,:,0] = c_root
    wake_network[1,:,0] = 2.0*c_root
    wake_network[:,:,1] = np.linspace(0.0, y_tip, N_spanwise)[np.newaxis,:]
    networks.append(wake_network)
    network_names.append("wake")

    return networks, network_names


def pan_air_pressures_to_csv():
    # Reads in the PAN AIR output data and writes the pressure distributions to a csv file

    # Read in output
    with open(pan_air_output_file, 'r') as panair_output:
        lines = panair_output.readlines()

    # Initialize data storage
    header = "jc,ip,x,y,z,wx,wy,wz,cp2ndu,cpisnu,lmachu,source,doublet"
    data = []

    # Get data
    within_network = False
    for line in lines:
        
        # Check for start of network
        if "network id" in line and "source type" in line:
            within_network = True
            continue

        # Check for end of network
        if "force / moment data for network" in line:
            within_network = False
            continue

        # If we're in a network, try to parse the data
        if within_network:
            split_line = line.split()
            if len(split_line) < 13:
                continue
            if split_line[0] != "jc":
                data.append(split_line)

    # Write out data
    with open(pan_air_data_file, 'w') as csv_handle:
        print(header, file=csv_handle)
        for data_line in data:
            print(",".join(data_line), file=csv_handle)


def create_vtk(N_chordwise_fore, N_chordwise_aft, N_spanwise, grid_label, split_dir="right"):
    # Creates delta wing mesh

    # Determine total number of vertices
    N_around_chord = (N_chordwise_fore + N_chordwise_aft - 1)*2 - 2
    N_verts = N_around_chord*(N_spanwise - 1) + 1
    vertices = np.zeros((N_verts,3))

    # Determine distributions of spanwise locations
    Y = np.linspace(0.0, b_semi, N_spanwise)
    c_of_y = np.linspace(c_root, 0.0, N_spanwise)
    x_le_of_y = np.linspace(0.0, c_root, N_spanwise)
    x_t_max_of_y = x_le_of_y + c_t_max*c_of_y
    z_t_max_of_y = np.linspace(0.5*t_max*c_root, 0.0, N_spanwise)

    # Loop through spanwise locations starting at the root to generate vertices
    for i, yi in enumerate(Y[:-1]):

        # Get current slice
        start = i*N_around_chord
        max_t_u = start + N_chordwise_aft
        le = max_t_u + N_chordwise_fore - 1
        max_t_l = le + N_chordwise_fore - 1
        end = (i+1)*N_around_chord

        # Set y coordinate
        vertices[start:end,1] = yi

        # Set x coordinate
        vertices[start:max_t_u,0] = np.linspace(c_root, x_t_max_of_y[i], N_chordwise_aft)
        vertices[max_t_u:le,0] = np.linspace(x_t_max_of_y[i], x_le_of_y[i], N_chordwise_fore)[1:]
        vertices[le:max_t_l,0] = np.linspace(x_le_of_y[i], x_t_max_of_y[i], N_chordwise_fore)[1:]
        vertices[max_t_l:end,0] = np.linspace(x_t_max_of_y[i], c_root, N_chordwise_aft)[1:-1]

        # Set z coordinate
        vertices[start:max_t_u,2] = np.linspace(0.0, z_t_max_of_y[i], N_chordwise_aft)
        vertices[max_t_u:le,2] = np.linspace(z_t_max_of_y[i], 0.0, N_chordwise_fore)[1:]
        vertices[le:max_t_l,2] = np.linspace(0.0, -z_t_max_of_y[i], N_chordwise_fore)[1:]
        vertices[max_t_l:end,2] = np.linspace(-z_t_max_of_y[i], 0.0, N_chordwise_aft)[1:-1]

    # Set tip
    vertices[-1,0] = c_root
    vertices[-1,1] = b_semi

    # Loop through spanwise locations starting at the root to generate panels
    panels = []
    for i in range(N_spanwise-2):

        # Determine start vertex index towards the root
        i_root_start = i*N_around_chord

        # Determine start vertex index towards the tip
        i_tip_start = (i+1)*N_around_chord

        # Loop around chord
        for j in range(N_around_chord):

            # Check we're not at the end
            if j != N_around_chord-1:

                if split_dir=="right": # The actual name of split_dir is arbitrary; I just had to name it something to distinguish

                    # Create first panel
                    panels.append([i_root_start+j, i_tip_start+j, i_root_start+j+1])

                    # Create second panel
                    panels.append([i_root_start+j+1, i_tip_start+j, i_tip_start+j+1])

                else:

                    # Create first panel
                    panels.append([i_root_start+j, i_tip_start+j, i_tip_start+j+1])

                    # Create second panel
                    panels.append([i_root_start+j+1, i_root_start+j, i_tip_start+j+1])

            # If we are at the end, loop back around
            else:

                if split_dir=="right":

                    # Create first panel
                    panels.append([i_root_start+j, i_tip_start+j, i_root_start])

                    # Create second panel
                    panels.append([i_root_start, i_tip_start+j, i_tip_start])

                else:

                    # Create first panel
                    panels.append([i_root_start+j, i_tip_start+j, i_tip_start])

                    # Create second panel
                    panels.append([i_root_start, i_root_start+j, i_tip_start])

    # Create final section at tip
    i_start = N_around_chord*(N_spanwise-2)
    for j in range(N_around_chord-1):
        panels.append([i_start+j, N_verts-1, i_start+j+1])
    panels.append([i_start+N_around_chord-1, N_verts-1, i_start])
    panels = np.array(panels)

    # Save
    mesh_file = study_dir + "meshes/delta_wing_{0}_{1}_split.vtk".format(grid_label, split_dir)
    _export_vtk(mesh_file, vertices, panels)

    return mesh_file


def run_machline(M, alpha, grid, mesh_file):
    # Runs MachLine at the given Mach number and angle of attack

    # Declare MachLine input
    case_name = mesh_file.replace(study_dir + "meshes/", "").replace(".vtk", "")
    body_file = study_dir + "results/" + case_name + "_{0}_{1}_aoa.vtk".format(grid, int(alpha))
    report_file = study_dir + "reports/" + case_name + "_{0}_{1}_aoa.json".format(grid, int(alpha))
    input_dict = {
        "flow": {
            "freestream_velocity": [np.cos(np.radians(alpha)), 0.0, np.sin(np.radians(alpha))],
            "gamma" : 1.4,
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "mirror_about" : "xz",
            "max_continuity_angle" : 1.0,
            "wake_model": {
                "wake_present" : True,
                "append_wake" : False
            },
            "reference": {
                "area": A_ref
            }
        },
        "solver": {
            "formulation": "morino",
            "matrix_solver": "GMRES",
            "run_checks": True
        },
        "post_processing" : {
            "pressure_rules" : {
                "second-order" : True,
                "isentropic" : True,
                "slender-body" : True,
                "linear" : True
            },
            "pressure_for_forces" : "second-order"
        },
        "output" : {
            "verbose" : False,
            "body_file" :  body_file,
            "report_file": report_file
        }
    }

    # Dump
    input_file = study_dir + "delta_wing_input.json"
    write_input_file(input_dict, input_file)

    # Run
    reports = run_quad(input_file, run=True)

    return reports


def compare(M, alpha, N_chordwise_fore, N_chordwise_aft, N_spanwise, grid, split_dir="right"):
    # Runs the delta wing through both MachLine and PAN AIR to compare

    # PAN AIR

    # Create networks
    networks, network_names = create_networks(N_chordwise_fore, N_chordwise_aft, N_spanwise)

    # Set up PAN AIR case
    case = pnr.PanairWrapper("Love Delta Wing",
                             description="Cory Goates, USU AeroLab",
                             directory="./"+study_dir)
    case._panair_loc = "./"+study_dir
    case.set_aero_state(mach=M, alpha=alpha, beta=0.0)

    # Add networks
    for network, name in zip(networks, network_names):
        if name=="wake":
            case.add_network(name, network, xy_indexing=False, network_type=18)
        else:
            case.add_network(name, network, xy_indexing=True)

    # Set extra info
    case.set_symmetry(True, False)
    case.set_reference_data(A_ref, b_semi, c_root)

    # Run PAN AIR
    result = case.run()

    # Extract forces
    FM = result.get_forces_and_moments()
    C_F_panair = np.zeros(3)
    C_F_panair[0] = FM["fx"]
    C_F_panair[1] = FM["fy"]
    C_F_panair[2] = FM["fz"]
    print("PAN AIR Forces: ", C_F_panair)

    # Extract data
    pan_air_pressures_to_csv()

    # MachLine

    # Create vtk file
    mesh_file = create_vtk(N_chordwise_fore, N_chordwise_aft, N_spanwise, grid, split_dir=split_dir)

    # Run MachLine
    reports = run_machline(M, alpha, grid, mesh_file)
    machline_forces = []
    for report, case in zip(reports, cases):
        try:
            machline_forces.append(report["total_forces"])
            print(case, report["total_forces"])
        except KeyError:
            machline_forces.append(None)
            continue

    # Plot pressures
    plot_pressure_slices(M, alpha, reports, grid, split_dir=split_dir)

    Cx = np.array([C_F_panair[0]] + [x["Cx"] for x in machline_forces])
    Cz = np.array([C_F_panair[2]] + [x["Cz"] for x in machline_forces])
    t = np.array([result.get_runtime()] + [x["total_runtime"] for x in reports])
    l_avg = reports[0]["mesh_info"]["average_characteristic_length"]

    return Cx, Cz, t, l_avg


def plot_pressure_slices(M, alpha, reports, grid, split_dir="right"):
    # Plots pressure slices from MachLine and PAN AIR

    # Get PAN AIR data
    header_pa, data_pa = get_data_from_csv(pan_air_data_file, remove_csv=False)
    x_pa = get_data_column_from_array(header_pa, data_pa, "x")
    y_pa = get_data_column_from_array(header_pa, data_pa, "y")
    C_P_pa = get_data_column_from_array(header_pa, data_pa, "cpisnu")

    # Get unique y stations out of PAN AIR data
    y_unique, ind_inv = np.unique(y_pa, return_inverse=True)

    # Loop through spanwise stations
    for i, y in enumerate(y_unique):

        percent_semi = int(round(y/b_semi, 2)*100)

        # Get indices of the data for this station
        ind = np.where(ind_inv == i)

        # Plot PAN AIR data
        plt.figure()
        plt.plot(x_pa[ind]/c_root, C_P_pa[ind], 'kv', label='PAN AIR', markersize=4)

        # Get data from MachLine and plot
        for report, case, line_style in zip(reports, cases, line_styles):

            # Get slice
            result_file = report["input"]["output"]["body_file"]
            headers_ml, data_ml = extract_plane_slice(result_file, [0.0, 1.0, 0.0], [0.0, y, 0.0], which_data='cell')
            x_ml = get_data_column_from_array(headers_ml, data_ml, "centroid:0")
            C_P_ml = get_data_column_from_array(headers_ml, data_ml, "C_p_ise")

            if alpha==0.0:

                # Figure out where the leading edge is so we only plot one surface
                i_le = len(x_ml)//2

                # If we're not at the tip (where PAN AIR has triangular panels), average the MachLine data to be apples-to-apples with PAN AIR
                if i < len(y_unique)-1:
                    x_ml_avg = 0.5*(x_ml[:i_le:2] + x_ml[1:i_le:2])
                    C_P_ml_avg = 0.5*(C_P_ml[:i_le:2] + C_P_ml[1:i_le:2])
                    plt.plot(x_ml_avg/c_root, C_P_ml_avg, line_style, label=case, linewidth=0.5)
                else:
                    plt.plot(x_ml[:i_le]/c_root, C_P_ml[:i_le], line_style, label=case, linewidth=0.5)

            else:

                # If we're not at the tip (where PAN AIR has triangular panels), average the MachLine data to be apples-to-apples with PAN AIR
                if i < len(y_unique)-1:
                    x_ml_avg = 0.5*(x_ml[::2] + x_ml[1::2])
                    C_P_ml_avg = 0.5*(C_P_ml[::2] + C_P_ml[1::2])
                    plt.plot(x_ml_avg/c_root, C_P_ml_avg, line_style, label=case, linewidth=0.5)
                else:
                    plt.plot(x_ml/c_root, C_P_ml, line_style, label=case, linewidth=0.5)

        # Format
        plt.xlabel('$x/c_r$')
        plt.ylabel('$C_P$')
        plt.gca().invert_yaxis()
        if alpha==0.0:
            plt.legend()
        plt.savefig(plot_dir + "pressure_at_alpha_{0}_percent_semispan_{1}_{2}_{3}_split.pdf".format(int(alpha), percent_semi, grid, split_dir))
        plt.savefig(plot_dir + "pressure_at_alpha_{0}_percent_semispan_{1}_{2}_{3}_split.svg".format(int(alpha), percent_semi, grid, split_dir))
        plt.close()


if __name__=="__main__":

    # Mesh parameters
    grids = ['coarse', 'medium', 'fine']
    Ncfs = [3, 6, 12]
    Ncas = [6, 12, 24]
    Nss = [10, 20, 40]
    #split_dirs = ["left", "right"]
    split_dirs = ["right"] # Verified on 6/22/23 that right is the best option here, as would be expected

    # Flow options
    alphas = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    M = 1.62

    # Loop
    for split_dir in split_dirs:

        # Initialize storage
        Cx = np.zeros((len(grids), len(alphas), 5))
        Cz = np.zeros((len(grids), len(alphas), 5))
        t = np.zeros((len(grids), len(alphas), 5))
        l_avg = np.zeros(len(grids))

        for i, grid in enumerate(grids):
            for j, alpha in enumerate(alphas):
                Cx[i,j,:], Cz[i,j,:], t[i,j,:], l_avg[i] = compare(M, alpha, Ncfs[i], Ncas[i], Nss[i], grid, split_dir=split_dir)

        # Plot axial force convergence
        plt.figure()
        for i, grid in enumerate(grids):
            if i==0:
                plt.plot(alphas, Cx[i,:,0], 'v', markersize=8-i*3, label="PAN AIR", mfc='none', mec='k')
            else:
                plt.plot(alphas, Cx[i,:,0], 'v', markersize=8-i*3, mfc='none', mec='k')
        plt.xlabel("$\\alpha$")
        plt.ylabel("$C_x$")
        plt.legend()
        plt.savefig(plot_dir + "C_x_convergence_{0}_split.pdf".format(split_dir))
        plt.savefig(plot_dir + "C_x_convergence_{0}_split.svg".format(split_dir))
        plt.close()

        # Plot normal force convergence
        plt.figure()
        for i, grid in enumerate(grids):
            if i==0:
                plt.plot(alphas, Cz[i,:,0], 'v', markersize=8-i*3, label="PAN AIR", mfc='none', mec='k')
            else:
                plt.plot(alphas, Cz[i,:,0], 'v', markersize=8-i*3, mfc='none', mec='k')
        plt.xlabel("$\\alpha$")
        plt.ylabel("$C_z$")
        plt.legend()
        plt.savefig(plot_dir + "C_z_convergence_{0}_split.pdf".format(split_dir))
        plt.savefig(plot_dir + "C_z_convergence_{0}_split.svg".format(split_dir))
        plt.close()

        for i, grid in enumerate(grids):

            # Plot axial force comparison
            plt.figure()
            plt.plot(alphas, Cx[i,:,0], 'kv', markersize=3, label="PAN AIR")
            for j, case in enumerate(cases):
                plt.plot(alphas, Cx[i,:,j+1], line_styles[j], label=case, linewidth=0.5)
            plt.xlabel("$\\alpha$")
            plt.ylabel("$C_x$")
            plt.legend()
            plt.savefig(plot_dir + "C_x_comparison_{0}_{1}.pdf".format(grid, split_dir))
            plt.savefig(plot_dir + "C_x_comparison_{0}_{1}.svg".format(grid, split_dir))
            plt.close()

            # Plot normal force comparison
            plt.figure()
            plt.plot(alphas, Cz[i,:,0], 'kv', markersize=3, label="PAN AIR")
            for j, case in enumerate(cases):
                plt.plot(alphas, Cz[i,:,j+1], line_styles[j], label=case, linewidth=0.5)
            plt.xlabel("$\\alpha$")
            plt.ylabel("$C_z$")
            plt.savefig(plot_dir + "C_z_comparison_{0}_{1}.pdf".format(grid, split_dir))
            plt.savefig(plot_dir + "C_z_comparison_{0}_{1}.svg".format(grid, split_dir))
            plt.close()

        # Calculate average run times over angles of attack
        t_avg = np.average(t, axis=1)
        t_std = np.std(t, axis=1)

        # Plot run times
        plt.figure()
        plt.plot(grids, t_avg[:,0], 'kv', markersize=3, label="PAN AIR")
        for j, case in enumerate(cases):
            plt.plot(grids, t_avg[:,j+1], line_styles[j], label=case, linewidth=0.5)
        plt.xlabel("Mesh Refinement")
        plt.ylabel("Total Runtime $[s]$")
        plt.yscale("log")
        plt.legend()
        plt.savefig(plot_dir + "time_comparison_{0}.pdf".format(split_dir))
        plt.savefig(plot_dir + "time_comparison_{0}.svg".format(split_dir))
        plt.close()

        # Calculate orders of convergence
        order_PA = np.zeros((2,len(alphas)))
        order_ML = np.zeros((2,len(alphas)))
        order_MH = np.zeros((2,len(alphas)))
        order_SL = np.zeros((2,len(alphas)))
        order_SH = np.zeros((2,len(alphas)))
        for i, alpha in enumerate(alphas):

            # Get PAN AIR orders of convergence
            order_PA[0,i] = get_order_of_convergence(l_avg*np.sqrt(2), Cx[:,i,0], truth_from_results=True)
            order_PA[1,i] = get_order_of_convergence(l_avg*np.sqrt(2), Cz[:,i,0], truth_from_results=True)

            # Get MachLine orders of convergence
            order_ML[0,i] = get_order_of_convergence(l_avg, Cx[:,i,1], truth_from_results=True)
            order_ML[1,i] = get_order_of_convergence(l_avg, Cz[:,i,1], truth_from_results=True)

            order_MH[0,i] = get_order_of_convergence(l_avg, Cx[:,i,2], truth_from_results=True)
            order_MH[1,i] = get_order_of_convergence(l_avg, Cz[:,i,2], truth_from_results=True)

            order_SL[0,i] = get_order_of_convergence(l_avg, Cx[:,i,3], truth_from_results=True)
            order_SL[1,i] = get_order_of_convergence(l_avg, Cz[:,i,3], truth_from_results=True)

            order_SH[0,i] = get_order_of_convergence(l_avg, Cx[:,i,4], truth_from_results=True)
            order_SH[1,i] = get_order_of_convergence(l_avg, Cz[:,i,4], truth_from_results=True)

        # Plot
        plt.figure()
        plt.plot(alphas, order_PA[0,:], 'kv', label="PAN AIR")
        plt.plot(alphas, order_ML[0,:], line_styles[0], label="ML")
        plt.plot(alphas, order_MH[0,:], line_styles[1], label="MH")
        plt.plot(alphas, order_SL[0,:], line_styles[2], label="SL")
        plt.plot(alphas, order_SH[0,:], line_styles[3], label="SH")
        plt.xlabel("$\\alpha$")
        plt.ylabel("Order of Convergence in $C_x$")
        plt.legend()
        plt.savefig(plot_dir + "Cx_convergence_order_{0}_split.pdf".format(split_dir))
        plt.savefig(plot_dir + "Cx_convergence_order_{0}_split.svg".format(split_dir))
        plt.close()

        plt.figure()
        plt.plot(alphas, order_PA[1,:], 'kv', label="PAN AIR")
        plt.plot(alphas, order_ML[1,:], line_styles[0], label="ML")
        plt.plot(alphas, order_MH[1,:], line_styles[1], label="MH")
        plt.plot(alphas, order_SL[1,:], line_styles[2], label="SL")
        plt.plot(alphas, order_SH[1,:], line_styles[3], label="SH")
        plt.xlabel("$\\alpha$")
        plt.ylabel("Order of Convergence in $C_z$")
        plt.savefig(plot_dir + "Cz_convergence_order_{0}_split.pdf".format(split_dir))
        plt.savefig(plot_dir + "Cz_convergence_order_{0}_split.svg".format(split_dir))
        plt.close()

        # Print out average orders of convergence
        print()
        print("Orders of Convergence")
        print("Case                   Cx                    Cz")
        print("-----------------------------------------------")
        print("PAN AIR", str(round(np.average(order_PA[0,:]), 3))+"+/-"+str(round(np.std(order_PA[0,:]), 3)), str(round(np.average(order_PA[1,:]), 3))+"+/-"+str(round(np.std(order_PA[1,:]), 3)))
        print("ML", str(round(np.average(order_ML[0,:]), 3))+"+/-"+str(round(np.std(order_ML[0,:]), 3)), str(round(np.average(order_ML[1,:]), 3))+"+/-"+str(round(np.std(order_ML[1,:]), 3)))
        print("MH", str(round(np.average(order_MH[0,:]), 3))+"+/-"+str(round(np.std(order_MH[0,:]), 3)), str(round(np.average(order_MH[1,:]), 3))+"+/-"+str(round(np.std(order_MH[1,:]), 3)))
        print("SL", str(round(np.average(order_SL[0,:]), 3))+"+/-"+str(round(np.std(order_SL[0,:]), 3)), str(round(np.average(order_SL[1,:]), 3))+"+/-"+str(round(np.std(order_SL[1,:]), 3)))
        print("SH", str(round(np.average(order_SH[0,:]), 3))+"+/-"+str(round(np.std(order_SH[0,:]), 3)), str(round(np.average(order_SH[1,:]), 3))+"+/-"+str(round(np.std(order_SH[1,:]), 3)))