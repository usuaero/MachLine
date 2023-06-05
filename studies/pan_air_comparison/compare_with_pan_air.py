import numpy as np
import panairwrapper.panairwrapper as pnr

from scipy.spatial import ConvexHull
from dev.helper_scripts.geometry_creator import _export_vtk

# Geometry parameters (for a Love's delta wing #11)
t_max = 0.08
c_t_max = 0.18
c_root = 1.0
b_semi = 1.0


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


def compare(M, alpha, N_chordwise_fore, N_chordwise_aft, N_spanwise):
    # Runs the delta wing through both MachLine and PAN AIR to compare

    # Create networks
    networks, network_names = create_networks(N_chordwise_fore, N_chordwise_aft, N_spanwise)

    # Set up PAN AIR case
    case = pnr.PanairWrapper("Love Delta Wing",
                             description="Cory Goates, USU AeroLab",
                             directory="./studies/pan_air_comparison/")
    case._panair_loc = "./studies/pan_air_comparison/"
    case.set_aero_state(mach=M, alpha=alpha, beta=0.0)

    # Add networks
    for network, name in zip(networks, network_names):
        case.add_network(name, network)

    # Set symmetry
    case.set_symmetry(True, True)

    # Run PAN AIR
    result = case.run()

    print(result)
    print(result.get_forces_and_moments())


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