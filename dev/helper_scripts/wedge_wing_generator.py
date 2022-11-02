# This script generates a vtk mesh file for symmetric wedge wings with an arbitrary taper ratio
# Airfoil coordinates are used throughout this script

from urllib.robotparser import RobotFileParser
import numpy as np
import matplotlib.pyplot as plt
import copy
import os
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
    Leading_edge_x = Zoc * np.tan(sweep_angle)
    TE_slope = (semispan * np.tan(sweep_angle) + taper - 1) / semispan
    TE_xloc = TE_slope * Zoc + 1

    return Zoc, Leading_edge_x, TE_xloc


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
    t_c = kwargs.get('t_c', 0.15) / 2 # The thickness is split in half here because 
                                      # it will be mirrored about the chordline for the lower surface
    node_ratio = kwargs.get('node_shift', 0.5)
    chordwise_cluster = kwargs.get('cosine_cluster', True)
    local_chord = kwargs.get('local_chord', 1.0)

    # Determine number of nodes, preventing needle panels near tips
    chord_nodes = int(np.ceil(local_chord / chord_spacing))
    # chord_nodes = cw_nodes
    if chord_nodes < 5:
        chord_nodes = 5
    # Check if at the wingtip where there should only be one node
    if abs(local_chord) < 10e-12:
        chord_nodes = 1
        x_vert = np.array([0.0], dtype=object)
        y_vert = np.array([0.0], dtype=object)
        return x_vert, y_vert
        
    # Identify ratio of nodes forward and aft of max thickness location
    LE_nodes = max(round(chord_nodes / (1/node_ratio)), 1)
    aft_nodes = max(chord_nodes - LE_nodes, 2)

    # Account for when we remove the duplicate value
    LE_nodes += 1

    # Cluster chorwise direction along the x axis
    if chordwise_cluster:
        # Clustering forward of max thickness location
        vTheta = np.linspace(0,np.pi, LE_nodes)
        fwd_xloc = 0.5 * (x_c) * (1-np.cos(vTheta))

        # Clustering aft of the max thickness location
        vTheta = np.linspace(0, np.pi, aft_nodes)
        aft_xloc = 0.5 * (1.0-x_c) * (1.0-np.cos(vTheta)) + x_c
        

    # No chordwise clustering
    else:
        fwd_xloc = np.linspace(0, x_c, LE_nodes)
        aft_xloc = np.linspace(x_c, 1.0, aft_nodes)
    

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
    LE_xloc = kwargs.get('LE_xloc', 0.0)
    local_chord = kwargs.get('local_chord', 1.0)
    
    # Scale airfoil and shift based on local chord length and leading edge x location
    x_vert = x_root * local_chord + LE_xloc 
    y_vert = y_root * local_chord

    breakpoint()
    print()
    

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
    plane : str
        plane in which to measure the 2D distance ex. 'xz', 'yz', 'xz'
    dimension : str
        select '2D' or '3D'

    Returns
    -------
    distance : float
        2-Dimensional distance between the two vertices along the airfoil plane
    """

    dim = kwargs.get('dimension', '3D')

    if dim == '2D':
        plane = kwargs.get('plane', 'xz')

        if plane == 'xy':
            dist = np.sqrt((v2[1]-v1[1])**2 + (v2[0]-v1[0])**2) 
        elif plane == 'xz':
            dist = np.sqrt((v2[2]-v1[2])**2  + (v2[0]-v1[0])**2) 
        elif plane == 'yz':
            dist = np.sqrt((v2[2]-v1[2])**2 + (v2[1]-v1[1])**2) 
        else:
            print(plane, ' is an unrecognized entry for the airfoil plane. Quitting')
            exit()
    else:
        dist = np.sqrt((v2[2] - v1[2])**2 + (v2[1]-v1[1])**2 + (v2[0]-v1[0])**2)

    return dist

def between(value, bound_1, bound_2):
    """Checks if a value is between two boundary values.
    Parameters
    ----------
    bound_1 : int
        boundary for check. Can be greater or less than boundary 2
    bound_2 : int
        boundary for check. Can be greater or less than boundary 1
    value : int
        value to verify if between bound_1 and bound_2
        
    Returns
    -------
    check : True/False
        boolean identifying if the passed value is between the 2 given boundaries
    """
    if (value == bound_1) and (value == bound_2):
        print("All values when comparing distances cannot be the same. Quitting...")
        exit()

    if bound_1 > bound_2:
        check = value in list(range(bound_2, bound_1+1))
    elif bound_2 > bound_1:
        check = value in list(range(bound_1, bound_2+1))
    else:
        check = False

    return check


def create_panels(airfoil_locations, vertices, le_list, xc_max):
    """Creates a list of panels, using ccw rrh to assign order of vertices.
    
    Parameters
    ----------
    airfoil_locations : list of lits 
        Contains the vertex locations of the TE, LE, and final vertex of each airfoil location
    vertices : list of lists
        Each inner list contains the x, y, z coordinates of the vertex locations
    le_list : list
        Single list containing the vertex locations of the leading edge vertices
    xc_max : list
        Single list containing the vertex locations of the max thickness location vertices
    
    Returns
    -------
    panels : list of lists
        Each inner list is a panel, identifying the orientation and location of each panel along the mesh
    """

    # Initialize panel container and number of vertices
    panels = []

    # Create a modified list of airfoil locations over upper surface
    upper_surf = copy.deepcopy(airfoil_locations)
    
    for semi in upper_surf:
        semi.remove(semi[-1])

    # Iterate over airfoil locations along upper surface
    for i, root in enumerate(upper_surf):
        if i < len(airfoil_locations)-1:
            tip = upper_surf[i+1]

            root_list = list(range(root[0], root[1]+1))
            tip_list = list(range(tip[0], tip[1]+1)) 


            # Create panels over upper surface of wing
            panel_iter = panel_generator(i, airfoil_locations, root_list, tip_list, vertices, le_list, xc_max, 'upper')
            panels += panel_iter


    # Create a modified list of airfoil locations over the lower surface
    lower_surf = copy.deepcopy(airfoil_locations)
    for semi in lower_surf:
        semi.remove(semi[0])


    for i, root in enumerate(lower_surf):
        if i < len(airfoil_locations)-1:
            tip = lower_surf[i+1]

            root_list = list(range(root[0], root[1]+1))
            tip_list = list(range(tip[0], tip[1]+1))

            # Create panels over lower surface of wing
            panel_iter = panel_generator(i, airfoil_locations, root_list, tip_list, vertices, le_list, xc_max, 'lower')
            panels += panel_iter

            # Connect final panel on lower surface to trailing edge
            p_one = airfoil_locations[i][0]
            p_two = airfoil_locations[i+1][0]

            # Determine the middle vertex based on what isn't already used
            if p_one not in panels[-1]:
                p_three = airfoil_locations[i][-1]
            else:
                p_three = airfoil_locations[i+1][-1]

            # Avoid zero area panels
            if not ((p_two == p_three) or (p_three == p_one)):
            
                panels.append([p_one, p_two, p_three]) # one, two, three

    return panels


def panel_generator(iter, airfoil_locations, root, tip, vertices, le_list, xc_max, upper_or_lower):

    # Initialize counters and max iteration trackers
    i_root_max = int(abs(root[-1]-root[0])) + (upper_or_lower == 'lower') # Lower surface iterations requires 1 additional tip iteration
    i_tip_max = int(abs(tip[-1] - tip[0])) + 1
    i_root = 0
    i_tip = 0
    panel_iter = []


    # Iterate over vertices
    while (i_root < i_root_max) and (i_tip < i_tip_max): # Needs to be an and statement

        # Always start panel numbering with given vertex on next semispan location
        p_three = tip[i_tip]
        p_one = root[i_root]
        p_two = -1 # Default value to ensure it has been assigned properly
        # if (p_one ==54) or (p_three==54):
        #     breakpoint()
        # Run a series of checks to ensure proper choice of which vertex will be chosen for this iteration
        # Verify that the next tip and root locations are not on the next airfoil over
        if ((i_tip + 1) >= len(tip)):

            # Check if the next root location is on trailing edge
            if i_root <= (len(root)-2):
                
                next_root = root[i_root + 1]
                next_tip = airfoil_locations[iter+1][0]
            else:
                next_root = airfoil_locations[iter][0]
                next_tip = airfoil_locations[iter+1][0]
            

        elif ((i_root + 1) >= len(root)):
            next_root = airfoil_locations[iter][0]
            next_tip = tip[i_tip + 1]
        

        # Set default for  next_root and next_tip
        else:
            next_tip = tip[i_tip + 1]
            next_root = root[i_root + 1]


        # Check to see if next_tip will be the only node at the trailing edge
        if (next_tip in le_list) and (next_tip in xc_max):
            p_two = next_root
            i_root += 1
            
        # Check to see if next vertex is on leading edge
        elif (next_root in le_list) and (p_three in le_list):
            p_two = next_root
            i_root += 1


        elif (next_tip in le_list) and (p_one in le_list):
            p_two = next_tip
            i_tip += 1


        # Check to see if a vertex is on max thickness ridge
        elif (p_one in xc_max):
            # Check relation of tip vertex relative to next_root
            if upper_or_lower == 'lower':
                crit_location = le_list[iter+1]
            else:
                crit_location = tip[0]

            tip_max_xc = xc_max[xc_max.index(p_one) + 2]

            if crit_location == next_tip:
                pass
            # If next tip has passed max thickness location
            elif between(next_tip, crit_location, tip_max_xc):
                p_two = next_tip
                i_tip += 1           

           
        elif (p_three in xc_max):
            # Check relation of root vertex relative to next_tip
            if upper_or_lower == 'lower':
                crit_location = le_list[iter]
            else:
                crit_location = root[0]

            root_max_xc = xc_max[xc_max.index(p_three) - 2]

            # Account for case when comparing leading edge to trailing edge vertices as next tip point
            if next_tip == crit_location:
                pass

            # If next root has passed max thickness location
            elif between(next_root, crit_location, root_max_xc):
                p_two = next_root
                i_root += 1


        if p_two == -1:
            # Check for tip conidition with one vertex and measure distances
            dist_root = distance(vertices[p_three], vertices[next_root])
            dist_tip = distance(vertices[p_one], vertices[next_tip])            

            # Check if shortest distance is along tip airfoil
            if dist_tip < dist_root:
                # Increment along tip airfoil
                p_two = next_tip
                i_tip += 1
                
            else:
                # Increment along root airfoil
                p_two = next_root
                i_root += 1

        # Check if next point is the wingtip with a single point
        if (p_three == p_two) or (p_two == p_one):
            return panel_iter
        else:   
            # Append panels to list in ccw direction
            panel_iter.append([p_one, p_two, p_three]) # one, two, three

    return panel_iter


def wingtip_cap(vertices, airfoil_locations, le_locations, max_thickness, **kwargs):
    """ Creates caps for the wingtips on both sides of the mesh
    
    Parameters
    ----------
    cap_type : str
        type of cap to add to wingtip: edge
    scale_factor : float
        percent of semispan length that the wingtip cap protrudes
    cap_angle : float
        degree angle of cap leading edge relative to wing leading edge

    Returns
    -------
    panels : list of lists
        updated panels list that includes the wingtip cap panels
    vertices : list of lists
        updated vertex list that inclues the wingtip vertices
    le_locations : list
        updated list of leading edge locations
    max_thickness : list
        updated list of max thickness locations
    """
    # Pull in given key word arguments
    cap_type = kwargs.get('cap_type', 'flat')
    scale_factor = kwargs.get('scale_factor', 0.02)
    cap_angle = kwargs.get('cap_angle', 60.0) * np.pi / 180.0

    # Initialize containers
    cap_panels = []
    tip = airfoil_locations[-1]

    # Determine x location of each vertex along the cap
    nodes = tip[-1] - tip[0] + 1

    # Determine how far out cap will protrude in the z direction
    # if cap_type == 'edge':
    #     z_cap = vertices[airfoil_locations[-1][-1]][-1] * (1 + scale_factor)

    if cap_type == 'flat':
        z_cap = vertices[tip[-1]][-1]



        # Iterate over nodes at wingtip
        temp = tip[0]
        tip_verts = vertices[temp::]
        # for i, val in enumerate(tip_verts):
        upper = list(range(tip[0], tip[1]+1))
        lower = list(range(tip[2], tip[1],-1))

        cap_bounds = [[tip[0], tip[1]], [tip[1], tip[2]]] # Serves as airfoil locations for the cap
        

    else:
        print('Only "flat" cap types are supported at this time. Quitting...')
        exit()


    # Create dummy variables to avoid errors when generating panels
    le_list = [5,6,7]
    xc_max = [0,1,2]
    iter = 0
    
    # Create panels over the cap of the wingtip
    panel_iter = panel_generator(iter, cap_bounds, upper, lower, vertices, le_list, xc_max, 'upper')
    cap_panels += panel_iter


    return cap_panels

def mirror_body(vertices, panels):
    """Mirror vertices about the xy plane
    
    Parameters
    ----------
    vertices : list of lists
        Complete list of vertices used for non-mirrored body
    panels : list of lists
        Complete list of panles for non-mirrored body.    

    Returns
    -------
    combined_vertices : list of lists
        Complete list of vertices, mirrored and original body
    combined_panels : 2-D numpy array
        Complete list of panels, mirrored and original body
    """
    # Initialize containers
    combined_vertices = []
    combined_panels = []
    N = len(vertices)


    # Mirror vertices accross the xy plane
    mirrored_vertices = []
    for val in vertices:
        mirrored_vertices.append([val[0], val[1], -val[2]])

    # Mirror panels, correcting ordering for outward normal and referencing mirrored vertices
    mirrored_panels = []
    for pan in panels:
        mirrored_panels.append([pan[2]+N, pan[1]+N, pan[0]+N])

    # Combined mirrored panels and mirrored vertices
    combined_panels = np.array(np.ndarray.tolist(panels) + mirrored_panels, dtype=int)
    combined_vertices = vertices + mirrored_vertices

    return combined_vertices, combined_panels



if __name__ == '__main__':

    # Define necessary portions of wing geometry
    x_cr = 0.18 # nondimensional location of max thickness at root
    x_ct = 0.18 # nondimensional location of max thickness at tip
    tc = 0.08 # nondimensional max thickness
    b_2c = 0.2315 / 0.230 # semispan nondimensionalized by root chord
    R_T = 0.25 # taper ratio c_t/c_R
    LE_sweep = 44.85 # leading edge sweep angle in degrees
    node_ratio = .35 # User input ratio of nodes to be placed forward of the max thickness location
    mirror_xy = True # Mirrors body across xy plane

    # Initialize number of nodes in chord and spanwize directions
    cw_nodes = 30 # Number of nodes along the upper surface
    sw_nodes = 30 # Number of nodes along the semispan
    cluster_cw = True
    cluster_sw = False

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
    
    # breakpoint()
    # Iterate over semispan locations, updating the vertex information at each point
    for i, z in enumerate(zoc):
        # Store location of beginning of airfoil location
        airfoil_locs.append([vertex_cnt+1])

        c_local = te_xloc[i] - le_xloc[i]
        # c_local = 1.0
        x_coord, y_coord = init_airfoil(S_avg, local_chord=c_local, x_c=xc_local[i], t_c=tc, node_shift=node_ratio, cosine_cluster=cluster_cw)
        x_scaled, y_scaled = scale_airfoil(x_coord, y_coord, local_chord=c_local, LE_xloc=le_xloc[i])
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
    

    # Determine panels for non-mirrored mesh
    panels = create_panels(airfoil_locs, vertices, leading_edges, max_thickness)

    # Add caps to the wingtips
    if R_T != 0.0:
        cap_panels = wingtip_cap(vertices, airfoil_locs, leading_edges, max_thickness)
        
        panels += cap_panels
    

    # Mirror (if selected) and output final mesh to selected location
    if 'helper_scripts' in os.getcwd(): 
        os.chdir('../../')
    elif 'dev' in os.getcwd():
        os.chdir('../')

    # filename= 'studies/delta_wing/meshes/delta_wing_clustered_mesh.vtk'
    filename= 'studies/delta_wing/meshes/tapered_wing_test.vtk'


    panels = np.array(panels, dtype=int)

    if mirror_xy:
        combined_verts, combined_panels = mirror_body(vertices, panels)

        _export_vtk(filename, combined_verts, combined_panels)
    else:
        _export_vtk(filename, vertices, panels)

    print("Geometry created")
    print("file location: ", filename)    