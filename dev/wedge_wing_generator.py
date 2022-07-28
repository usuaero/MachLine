# This script generates a vtk mesh file for symmetric wedge wings with an arbitrary taper ratio
# Airfoil coordinates are used throughout this script

import numpy as np
import matplotlib.pyplot as plt


def airfoil_coord(x_c, t_c, node_ratio, z_location, chordwise_cluster=True):

    # Identify ratio of nodes forward of max thickness location
    if type(node_ratio) != float:
        LE_nodes = round(cw_nodes / (1/x_c))
    else:
        LE_nodes = round(cw_nodes / (1/node_ratio))
    
    aft_nodes = cw_nodes - LE_nodes
    # Account for when we remove the duplicate value
    LE_nodes += 1

    # Cluster chorwise direction along the x axis
    if chordwise_cluster:
        # Clustering forward of max thickness location
        vTheta = np.linspace(0,np.pi, LE_nodes)
        fwd_xloc = 0.5 * x_c * (1-np.cos(vTheta))

        # Clustering aft of the max thickness location
        vTheta = np.linspace(0, np.pi, aft_nodes)
        aft_xloc = 0.5 * (1-x_c) * (1-np.cos(vTheta)) + x_c

    # No chordwise clustering
    else:
        fwd_xloc = np.linspace(0, x_c, LE_nodes)
        aft_xloc = np.linspace(x_c, 1, aft_nodes)

    # Calculate the y location based on point slope form
    fwd_yloc = t_c / x_c * fwd_xloc # equation derived from point slope form
    aft_yloc = -t_c/(1-x_c) * aft_xloc + (t_c / (1 - x_c)) # equation derived from point slope form
    
    # Append fwd and aft wedge datapoints with removing redundant points
    fwd_xloc = fwd_xloc[:-1]
    fwd_yloc = fwd_yloc[:-1]
    xloc = np.append(fwd_xloc, aft_xloc)
    yloc = np.append(fwd_yloc, aft_yloc)

    return xloc, yloc

def spanwise_coord(spanwise_nodes, semispan, LE_sweep, R_T, spanwise_cluster=True):

    # Convert degrees to radians
    LE_sweep = LE_sweep * np.pi/180

    # Identify node locations along semispan
    # Cosine cluster points
    if spanwise_cluster:
        vTheta = np.linspace(0, np.pi, spanwise_nodes)
        Zob = 0.5 * semispan *  (1 - np.cos(vTheta))
    else:
        Zob = np.linspace(0, semispan, spanwise_nodes)

    # Determine x location at each z

    LE_xloc = Zob * np.tan(LE_sweep)
    TE_xloc = (semispan * np.tan(LE_sweep) + R_T - 1) / semispan * (Zob - 1)

    breakpoint()

    print()

    return Zob, LE_xloc, TE_xloc

def plot_wing_3D(x, y, z):

    ax = plt.axes(projection='3d')
    ax.scatter3D(x, y, z)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.view_init(30, 0) # pass angles to view the plot in degrees

    plt.show()

    # https://www.youtube.com/watch?v=gqoLLGgbeAE for seeing tutorial on plotting in 3D, including a wireframe

if __name__ == '__main__':

    # Define necessary portions of wing geometry
    xc = 0.18 # nondimensional location of max thickness
    tc = 0.08 # nondimensional max thickness
    b_2c = 1.0065 # semispan nondimensionalized by root chord
    R_T = 0.25 # taper ratio c_t/c_R
    LE_sweep = 44.85 # leading edge sweep angle in degrees
    node_ratio = "standard" # User input ratio of nodes to be placed forward of the max thickness location
                            # Default is number of nodes * xc rounded.


    # Initialize number of nodes in chord and spanwize directions
    cw_nodes = 40
    sw_nodes = 40
    cluster_cw = True
    cluster_sw = True



    Zoc = 0
    # Determine airfoil coordinates at root
    zloc, yloc = airfoil_coord(xc, tc, node_ratio, Zoc, cluster_cw)

    spanwise_coord(sw_nodes, b_2c, R_T, cluster_sw)
    
    # Calculate location of tip leading edge
    tip_LE = b_2c * np.tan(LE_sweep)
    tip_TE = tip_LE + R_T
    
    
    # Arrange vertices in ccw direction starting at root leading edge
    vertices = np.array([[0, 1, tip_TE, tip_LE, 0], # x axis
                         [0, 0, 0, 0, 0], # y axis
                         [-0, -0, -b_2c, -b_2c, -0]]) # z axis

    # Check planform outline
    plt.figure()
    # plt.scatter(xloc, yloc, c='k') # plot wedge airfoil
    # plt.plot(xloc, yloc, color='k', marker='.')
    plt.plot(vertices[2], vertices[0], 'k') # planform plot
    plt.plot([0.0, min(vertices[2])], [xc, tip_LE], 'k') # wedge airfoil ridge line
    plt.gca().invert_yaxis()  
    plt.gca().set_aspect('equal')
    plt.show()
    
    # Cosine cluster in the spanwise direction
    if cluster_sw:
        vTheta = np.linspace(0, np.pi, sw_nodes)

    else:
        sw_loc = np.linspace(0, b_2c, sw_nodes)

    

    print("Geometry created")