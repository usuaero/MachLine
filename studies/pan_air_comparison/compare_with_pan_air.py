import numpy as np
import panairwrapper.panairwrapper as pnr

from scipy.spatial import ConvexHull
from dev.helper_scripts.geometry_creator import _export_vtk

# Geometry parameters (for a Love's delta wing #11)
t_max = 0.08
c_t_max = 0.18
c_root = 0.23
b_semi = 0.463/2.0
A_ref = c_root*b_semi

# Run parameters
study_dir = "studies/pan_air_comparison/"
mesh_file = study_dir + "meshes/delta_wing.vtk"
body_file = study_dir + "results/delta_wing.vtk"
report_file = study_dir + "reports/report.json"


def create_networks(N_chordwise_fore, N_chordwise_aft, N_spanwise):
    # Creates the networks of vertices for PAN AIR to use

    networks = []
    network_names = []

    # Initialize fore network
    fore_network = np.zeros((N_spanwise, N_chordwise_fore, 3))

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
        fore_network[:,i,0] = np.linspace(x_root[i], x_tip, N_spanwise)
        fore_network[:,i,1] = np.linspace(y_root[i], y_tip, N_spanwise)
        fore_network[:,i,2] = np.linspace(z_root[i], z_tip, N_spanwise)

    # Add
    networks.append(fore_network)
    network_names.append("fore")

    # Initialize aft network
    aft_network = np.zeros((N_spanwise, N_chordwise_aft, 3))

    # Get points at root
    x_root = np.linspace(c_t_max*c_root, c_root, N_chordwise_aft)
    y_root = np.zeros(N_chordwise_aft)
    z_root = np.linspace(0.5*t_max*c_root, 0.0, N_chordwise_aft)

    # Interpolate
    for i in range(N_chordwise_aft):
        aft_network[:,i,0] = np.linspace(x_root[i], x_tip, N_spanwise)
        aft_network[:,i,1] = np.linspace(y_root[i], y_tip, N_spanwise)
        aft_network[:,i,2] = np.linspace(z_root[i], z_tip, N_spanwise)

    # Add
    networks.append(aft_network)
    network_names.append("aft")

    return networks, network_names


def create_vtk(N_chordwise_fore, N_chordwise_aft, N_spanwise):
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

                # Create first panel
                panels.append([i_root_start+j, i_tip_start+j, i_root_start+j+1])

                # Create second panel
                panels.append([i_root_start+j+1, i_tip_start+j, i_tip_start+j+1])

            # If we are at the end, loop back around
            else:

                # Create first panel
                panels.append([i_root_start+j, i_tip_start+j, i_root_start])

                # Create second panel
                panels.append([i_root_start, i_tip_start+j, i_tip_start])

    # Create final section at tip
    i_start = N_around_chord*(N_spanwise-2)
    for j in range(N_around_chord-1):
        panels.append([i_start+j, N_verts-1, i_start+j+1])
    panels.append([i_start+N_around_chord-1, N_verts-1, i_start])
    panels = np.array(panels)

    # Save
    _export_vtk(mesh_file, vertices, panels)


def compare(M, alpha, N_chordwise_fore, N_chordwise_aft, N_spanwise):
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
        case.add_network(name, network)

    # Set extra info
    case.set_symmetry(True, True)
    case.set_reference_data(A_ref, b_semi, c_root)

    # Run PAN AIR
    result = case.run()

    # Extract forces
    FM = result.get_forces_and_moments()
    C_F_panair = np.zeros(3)
    C_F_panair[0] = FM["fx"]
    C_F_panair[1] = FM["fy"]
    C_F_panair[2] = FM["fz"]

    # MachLine

    # Create vtk file
    create_vtk(N_chordwise_fore, N_chordwise_aft, N_spanwise)


if __name__=="__main__":

    # Mesh parameters
    grids = ['coarse']#, 'medium', 'fine']
    Ncfs = [5, 10, 20]
    Ncas = [10, 20, 40]
    Nss = [10, 20, 40]

    # Angles of attack
    alphas = [0.0]
    M = 1.62

    # Loop
    for i, grid in enumerate(grids):
        for j, alpha in enumerate(alphas):
            compare(M, alpha, Ncfs[i], Ncas[i], Nss[i])