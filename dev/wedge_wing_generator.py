# This script generates a vtk mesh file for symmetric wedge wings with an arbitrary taper ratio
# Airfoil coordinates are used throughout this script

import numpy as np
import matplotlib.pyplot as plt

def spanwise_coord(spanwise_nodes, semispan, **kwargs):
    """
    Determines z coordinates of nodes along the wingspan.
    Also determines the leading and trailing edges of the wing.

    Parameters
    ----------
    sweep_angle : float
        Leading edge sweep angle, given in radians
    taper_ratio : float
        Taper ratio
    cosine_cluster: bool
        Cosine clusters the spanwise nodes

    Output
    ------
    Zoc : numpy 1-D array
        Spanwise z locations of the nodes. Nondimensionalized by root chord
    LE_xloc : numpy 1-D array
        x location coordinate points of the leading edge
    TE_xloc : numpy 1-D array
        y location coordinate points of the trailing edge
    """

    sweep_angle = kwargs.get('sweep_angle', 0.0)
    taper = kwargs.get('taper_ratio', 1.0)
    cosine_cluster = kwargs.get('cosine_cluster', True)

    # Identify node locations along semispan
    # Cosine cluster points
    if cosine_cluster:
        vTheta = np.linspace(0, np.pi, spanwise_nodes)
        Zoc = 0.5 * semispan *  (1 - np.cos(vTheta))
    else:
        Zoc = np.linspace(0, semispan, spanwise_nodes)

    # Determine LE and TE locations at each z
    tip_LE = b_2c * np.tan(sweep_angle)
    tip_TE = tip_LE + R_T

    LE_xloc = Zoc * np.tan(sweep_angle)
    TE_slope = (semispan * np.tan(sweep_angle) + taper - 1) / semispan
    TE_xloc = TE_slope * Zoc + 1



    return Zoc, LE_xloc, TE_xloc


def init_airfoil(chord_spacing, **kwargs):
    """
    Returns the mesh coordinates along the xy plane at the wing root. (Airfoil coordinates)
    
    Parameters
    ----------
    x_c : float
        x location of max thickness(nondimensionalized by local chord length)
    t_c : float
        max thickness of airfoil (nondimensionalized by local chord length)
    LE_loc : float
        x location of leading edge of the airfoil (nondimensionalized by root chord)
    TE_loc : float
        x location of trailing edge of the airfoil (nondimensionalized by root chord)
    node_shift : float
        Ratio of nodes in front of max thickness location. Remainder will be aft of max thickness lcocation
    cosine_cluster : bool
        Determines if cosine clustering will be used for the chordwise distribution of node points

    Returns
    -------
    xloc : 1-D numpy array
        x locations associated with airfoil coordinates
    yloc : 1-D numpy array
        y locations associated with airfoil coordinates
    """
    # Pull in values from kwargs dict
    x_c = kwargs.get('x_c', 0.5)
    t_c = kwargs.get('t_c', 0.15)
    node_ratio = kwargs.get('node_shift', 0.5)
    chordwise_cluster = kwargs.get('cosine_cluster', True)

    # Determine number of nodes, preventing needle panels near tips
    chord_nodes = int(np.ceil(1 / chord_spacing))
    if chord_nodes < 5:
        chord_nodes += 1
        
    # Identify ratio of nodes forward and aft of max thickness location
    LE_nodes = round(chord_nodes / (1/node_ratio))
    aft_nodes = chord_nodes - LE_nodes

    # Account for when we remove the duplicate value
    LE_nodes += 1

    # Cluster chorwise direction along the x axis
    if chordwise_cluster:
        # Clustering forward of max thickness location
        vTheta = np.linspace(0,np.pi, LE_nodes)
        fwd_xloc = 0.5 * (x_c) * (1-np.cos(vTheta))

        # Clustering aft of the max thickness location
        vTheta = np.linspace(0, np.pi, aft_nodes)
        aft_xloc = 0.5 * (1-x_c) * (1-np.cos(vTheta)) + x_c

    # No chordwise clustering
    else:
        fwd_xloc = np.linspace(0, x_c, LE_nodes)
        aft_xloc = np.linspace(x_c, 1, aft_nodes)
    

    # Calculate the y location based on point slope form
    fwd_yloc = t_c / x_c * (fwd_xloc) # equation derived from point slope form
    aft_yloc = -t_c/(1-x_c) * (aft_xloc) + (t_c / (1 - x_c)) # equation derived from point slope form
    
    # Append fwd and aft wedge datapoints with removing redundant points
    fwd_xloc = fwd_xloc[:-1]
    fwd_yloc = fwd_yloc[:-1]
    x_vert = np.append(fwd_xloc, aft_xloc)
    y_vert = np.append(fwd_yloc, aft_yloc)

    return x_vert, y_vert


def scale_airfoil(x_root, y_root, **kwargs):
    """
    Returns the mesh coordinates along the xy plane at the wing root. (Airfoil coordinates)
    
    Parameters
    ----------
    x_root : numpy 1-D array
        x values for airfoil coordinates at root
    y_root : numpy 1-D array
        y values for airfoil coordinates at root
    x_c : float
        x location of max thickness(nondimensionalized by local chord length)
    t_c : float
        max thickness of airfoil (nondimensionalized by local chord length)
    LE_xloc : float
        x location of leading edge of the airfoil (nondimensionalized by root chord)
    TE_xloc : float
        x location of trailing edge of the airfoil (nondimensionalized by root chord)
    
    Returns
    -------
    x vertices: numpy 1-D array
        x vertices for scaled airfoil

    y vertices: numpy 1-D array
        y vertices for scaled airfoil
    """
    # Pull values from kwargs dictionary
    x_c = kwargs.get('x_c', 0.5)
    t_c = kwargs.get('t_c', 0.15)
    LE_xloc = kwargs.get('LE_xloc', 0.0)
    TE_xloc = kwargs.get('TE_xloc', 1.0)

    # Determine local chord relative to root chord
    c_local = TE_xloc - LE_xloc
    
    # Scale airfoil based on local chord size
    x_vert = x_root * c_local + LE_xloc # also shifts airfoil back to line up with local leading edge
    y_vert = y_root * c_local
    print(' Local chord:', c_local)
    return x_vert, y_vert

def plot_wing_3D(x, y, z):

    ax = plt.axes(projection='3d')
    for i,x_loc in enumerate(x):
        ax.scatter3D(x_loc, y[i], z)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    # ax.view_init(30, 0) # pass angles to view the plot in degrees

    

    # https://www.youtube.com/watch?v=gqoLLGgbeAE for seeing tutorial on plotting in 3D, including a wireframe

def plot_airfoil(x,y):
    fig, ax = plt.subplots()
    ax.scatter(x,y)
    ax.axis('equal')
    plt.show()

if __name__ == '__main__':

    # Define necessary portions of wing geometry
    x_cr = 0.18 # nondimensional location of max thickness at root
    x_ct = 0.18 # nondimensional location of max thickness at tip
    tc = 0.08 # nondimensional max thickness
    b_2c = 1.0065 # semispan nondimensionalized by root chord
    R_T = 0.0 # taper ratio c_t/c_R
    LE_sweep = 44.85 # leading edge sweep angle in degrees
    node_ratio = .25 # User input ratio of nodes to be placed forward of the max thickness location


    # Initialize number of nodes in chord and spanwize directions
    cw_nodes = 40
    sw_nodes = 20
    cluster_cw = True
    cluster_sw = True

    # Determine average node spacing at root chord
    S_avg = 1 / cw_nodes
    # Convert sweep angle to radians
    LE_sweep = LE_sweep * np.pi / 180

    # Determine spanwise locations where mesh will be generated
    zoc, le_xloc, te_xloc = spanwise_coord(sw_nodes, b_2c, sweep_angle=LE_sweep, taper_ratio=R_T, cosine_cluster=cluster_sw)
   
    # Interpolate the location of max thickness along the semispan of the wing
    xc_local = np.linspace(x_cr, x_ct, len(zoc))   
    
    # Determine airfoil coordinates at each root location
    x_coord, y_coord = init_airfoil(S_avg, x_c=xc_local[0], t_c=tc, node_shift=node_ratio, cosine_cluster=cluster_cw)
    plot_airfoil(x_coord, y_coord)

    for i, val in enumerate(zoc):
        print("Node ", i, end='')
        x_scaled, y_scaled = scale_airfoil(x_coord, y_coord, x_c=xc_local[i], t_c=tc, LE_xloc=le_xloc[i], TE_xloc=te_xloc[i])
        plot_airfoil(x_scaled, y_scaled)

        # plot_wing_3D(x_scaled, y_scaled, val)
    plt.show()

    
    # Calculate location of tip leading edge
    tip_LE = b_2c * np.tan(LE_sweep)
    tip_TE = tip_LE + R_T
    # breakpoint()
    # Plot airfoil at each semispan location
    # for i, semispan in enumerate(z_vert[i]):
        # plot_airfoil(x_vert[i], y_vert[i])

    # Arrange vertices in ccw direction starting at root leading edge
    vertices = np.array([[0, 1, tip_TE, tip_LE, 0], # x axis
                         [0, 0, 0, 0, 0], # y axis
                         [0, 0, b_2c, b_2c, 0]]) # z axis

    # Check planform outline
    plt.figure()
    plt.plot(vertices[2], vertices[0], 'k') # planform plot
    plt.plot([0.0, max(vertices[2])], [x_cr, tip_LE+x_ct*(tip_TE-tip_LE)], 'k') # wedge airfoil ridge line
    plt.gca().invert_yaxis()  # aligns x axis with airfoil coordinate system
    plt.gca().invert_xaxis()  # aligns z axis with airfoil coordinate system
    plt.gca().set_aspect('equal') # ensure plot isn't skewed in any direction
    plt.show()

    print("Geometry created")