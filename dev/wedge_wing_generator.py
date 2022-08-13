# This script generates a vtk mesh file for symmetric wedge wings with an arbitrary taper ratio
# Airfoil coordinates are used throughout this script

import numpy as np
import matplotlib.pyplot as plt
import copy
from sys import exit
from geometry_creator import _export_vtk


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
    # tip_LE = b_2c * np.tan(sweep_angle)
    # tip_TE = tip_LE + R_T

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
    local_chord : float
        local chord length (nondimensionalized by root chord length)
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
    local_chord = kwargs.get('local_chord', 1.0)

    # Determine number of nodes, preventing needle panels near tips
    chord_nodes = int(np.ceil(local_chord / chord_spacing))
    if chord_nodes < 3:
        chord_nodes = 3
    # Check if at the wingtip where there should only be one node
    if abs(local_chord) < 10e-12:
        chord_nodes = 1
        x_vert = np.array([0.0], dtype=object)
        y_vert = np.array([0.0], dtype=object)
        return x_vert, y_vert
        
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
    x_vert = np.append(aft_xloc[::-1], fwd_xloc[::-1]) # sort the data starting at traling edge
    y_vert = np.append(aft_yloc[::-1], fwd_yloc[::-1])

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
        leading edge x location of wing at identified semispan locations
    local_chord : float
        local chord length (nondimensionalized by root chord length)
    
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
    local_chord = kwargs.get('local_chord', 1.0)
    
    # Scale airfoil and shift based on local chord length and leading edge x location
    x_vert = x_root * local_chord + LE_xloc 
    y_vert = y_root * local_chord

    return x_vert, y_vert


def mirror_x(temp_list):
    """Mirrors a list without negating any results"""
    # Check if there is only 1 node
    if len(temp_list) == 1:
        return temp_list
    else:
        mirrored_list = np.append(temp_list[:-1], temp_list[::-1]) # remove duplicate leading edge vertex

        return mirrored_list[:-1] # remove duplicate trailing edge vertex

def mirror_y(temp_list):
    """Mirrors a list, negating the mirrored resutls"""
    # Check if there is only 1 node
    if len(temp_list) == 1:
        return temp_list
    else:

        mirrored_list = np.append(temp_list[:-1], -1*(temp_list[::-1])) # remove duplicate leading edge vertex

        return mirrored_list[:-1] # remove duplicate trailing edge vertex


def distance(v1, v2, **kwargs):
    """Calculates distance between two points, in either 2D or 3D.
    
    Parameters
    ----------
    airfoil_plane : str
        Plane which the airfoil sections lie along. ex 'xz', 'yz', 'xz'
    dimension : str
        select '2D' or '3D'

    Returns
    -------
    distance : float
        2-Dimensional distance between the two vertices along the airfoil plane
    """

    dim = kwargs.get('dimension', '3D')

    if dim == '2D':
        plane = kwargs.get('airfoil_plane', 'xy')

        if plane == 'xz':
            dist = np.sqrt((v2[1]-v1[1])**2 + (v2[0]-v1[0])**2) # compare xy distance
        elif plane == 'yz':
            dist = np.sqrt((v2[2]-v1[2])**2  + (v2[0]-v1[0])**2) # compare xz distance
        elif plane == 'xy':
            dist = np.sqrt((v2[2]-v1[2])**2 + (v2[1]-v1[1])**2) # compare yz distance
        else:
            print(plane, ' is an unrecognized entry for the airfoil plane. Quitting')
            exit()
    else:
        dist = np.sqrt((v2[2] - v1[2])**2 + (v2[1]-v1[1])**2 + (v2[0]-v1[0])**2)

    return dist


def create_panels(airfoil_locations, vertices, le_list, xc_max):
    """Creates a list of panels, using ccw rrh to assign order of vertices."""
    # Initialize panel container
    panels = []

    # Create a modified list of airfoil locations over upper surface
    upper_surf = copy.deepcopy(airfoil_locations)
    
    for semi in upper_surf:
        semi.remove(semi[-1])

    # Iterate over airfoil locations along upper surface
    for i, root in enumerate(upper_surf):
        if i < len(airfoil_locations)-1:
            tip = upper_surf[i+1]

            # Create panels over upper surface of wing
            panel_iter = panel_generator(root, tip, vertices, upper_surf, le_list, xc_max)
            panels += panel_iter


    # Create a modified list of airfoil locations over the lower surface
    lower_surf = copy.deepcopy(airfoil_locations)
    for semi in lower_surf:
        semi.remove(semi[0])


    for i, root in enumerate(lower_surf):
        if i < len(airfoil_locations)-1:
            tip = lower_surf[i+1]

            # Create panels over lower surface of wing
            panel_iter = panel_generator(root, tip, vertices, lower_surf, le_list, xc_max)
            panels += panel_iter

            # # Connect final panel on lower surface to trailing edge
            p_one = root[1] + 1
            p_three = root[1]

            # Determine the middle vertex based on what isn't already used
            if airfoil_locations[i][0] not in panels[-1]:
                p_two = airfoil_locations[i][0]
            else:
                p_two = airfoil_locations[i+1][0]
            
            panels.append([p_one, p_two, p_three])

                
    return panels


def panel_generator(root, tip, vertices, surface_verts, le_list, xc_max):

    # Iterating fwd over airfoil locations requires to exclude the last location
        
    # Initialize counters and max iteration trackers
    i_root_max = (root[1]-root[0]) + 1
    i_tip_max = (tip[1] - tip[0]) + 1
    i_root = 0
    i_tip = 0
    panel_iter = []
    error_iter = 0

    # Iterate over vertices
    while (i_root < i_root_max) and (i_tip < i_tip_max): # Needs to be an and statement

        # Always start panel numbering with given vertex on next semispan location
        p_one = tip[0] + i_tip
        p_three = root[0] + i_root



        # Minimize aspect ratio of panel
        next_tip = p_one + 1
        next_root = p_three + 1

        # Check to see if next_tip will be the only node at the trailing edge
        if (next_tip in le_list) and (next_tip in xc_max):
            p_two = p_three + 1
            i_root += 1
            
        # Check to see if next vertex is on leading edge
        elif (next_root in le_list) and ((p_one in le_list) or (p_three in le_list)):
            p_two = next_root
            i_tip += 1


        elif (next_tip in le_list) and ((p_one in le_list) or (p_three in le_list)):
            p_two = next_tip
            i_root += 1


        # Check to see if next vertex is on max thickness ridge
        elif (next_root in xc_max) and ((p_one in xc_max) or (p_three in xc_max)):
            p_two = next_root
            i_root += 1

            
        elif (next_tip in xc_max) and ((p_one in xc_max) or (p_three in xc_max)):
            p_two = next_tip
            i_tip += 1


        else:
            # Check for tip conidition with one vertex and measure distances
            max_vert = (surface_verts[-1][1]-1)
            if (p_one < max_vert) and (p_three < max_vert):
                # Determine diagonal distances
                # print('p_one: ', p_one)
                # print('p_three: ', p_three)
                # print('next_tip: ', next_tip)
                # print('next_root: ', next_root)
                dist_root = distance(vertices[p_one], vertices[next_root], dimension='2D')
                dist_tip = distance(vertices[p_three], vertices[next_tip], dimension='2D')

            # if only one node on tip airfoil
            else:
                dist_root = 0.0
                dist_tip = 1.0
                error_iter += 1
                if error_iter > 10:
                    print("You entered the never ending loop...quitting")
                    exit()
                # breakpoint()
            

            # Check if shortest distance is along tip airfoil
            if dist_tip < dist_root:
                # Increment along tip airfoil
                i_tip += 1
                p_two = tip[0] + i_tip
                
            else:
                # Increment along root airfoil
                i_root += 1
                p_two = root[0] + i_root

        if (p_one == p_two) or (p_two == p_three):
            return panel_iter
            

        # Append panels to list in ccw direction
        panel_iter.append([p_one, p_two, p_three])

    return panel_iter



def plot_airfoil(x,y):
    fig, ax = plt.subplots()
    ax.scatter(x,y)
    ax.set_xlim([0,1.])
    ax.set_ylim([-0.25, 0.25])
    ax.set_aspect('equal')
    plt.show()

def plot_planform(semispan, sweep_angle, taper, xcr, xct):
    # Calculate location of tip leading edge
    tip_LE = semispan * np.tan(sweep_angle)
    tip_TE = tip_LE + taper


    # Arrange vertices in ccw direction starting at root leading edge
    vertices = np.array([[0, 1, tip_TE, tip_LE, 0], # x axis
                         [0, 0, 0, 0, 0], # y axis
                         [0, 0, semispan, semispan, 0]]) # z axis

    # Check planform outline
    plt.figure()
    plt.plot(vertices[2], vertices[0], 'k') # planform plot
    plt.plot([0.0, max(vertices[2])], [xcr, tip_LE+xct*(tip_TE-tip_LE)], 'k') # wedge airfoil ridge line
    plt.gca().invert_yaxis()  # aligns x axis with airfoil coordinate system
    plt.gca().invert_xaxis()  # aligns z axis with airfoil coordinate system
    plt.gca().set_aspect('equal') # ensure plot isn't skewed in any direction
    plt.show()
    return

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
    cw_nodes = 20
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

    # Initialize vertex containers
    vertices = []
    airfoil_locs = []
    max_thickness = []
    leading_edges = []
    vertex_cnt = -1
    # Iterate over semispan locations, updating the vertex information at each point
   
    for i, z in enumerate(zoc):
        # Store location of beginning of airfoil location
        airfoil_locs.append([vertex_cnt+1])

        c_local = te_xloc[i] - le_xloc[i]
        x_coord, y_coord = init_airfoil(S_avg, local_chord=c_local, x_c=xc_local[0], t_c=tc, node_shift=node_ratio, cosine_cluster=cluster_cw)
        x_scaled, y_scaled = scale_airfoil(x_coord, y_coord, x_c=xc_local[i], t_c=tc, local_chord=c_local, LE_xloc=le_xloc[i])
        # mirror x and y scaled airfoils
        x_vert = mirror_x(x_scaled)
        y_vert = mirror_y(y_scaled)

        # Iterate over airfoil datapoints
        for j, x in enumerate(x_vert):
            vertices.append([x, y_vert[j], z])
            vertex_cnt += 1

            temp = (len(x_vert)//2 + len(x_vert)%2)
            # Store location of leading edge
            if (j == temp) or (temp==1):
                airfoil_locs[i].append(vertex_cnt)
            # Store max thickness location
            temp_thickness = c_local * xc_local[i] + le_xloc[i]
            if x == temp_thickness:
                max_thickness.append(vertex_cnt)
            if x == le_xloc[i]:
                leading_edges.append(vertex_cnt)
            
        
        # Store end location of airfoil location
        airfoil_locs[i].append(vertex_cnt)
    
    # Check vertices
    print("vertex count: ", vertex_cnt)
    max_dim = 0
    for vertex in vertices:
        if len(vertex) > max_dim:
            max_dim = len(vertex)
    
    print("Max dimension is ", max_dim)
    
    panels = np.array(create_panels(airfoil_locs, vertices, leading_edges, max_thickness), dtype=int)

    _export_vtk('wedge_wing_test.vtk', vertices, panels)




    print("Geometry created")