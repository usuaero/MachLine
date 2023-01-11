import numpy as np
from sys import exit
import pyvista as pv


def verify_xz_symmetry(points):

    # Iterate one point at a time
    for i, pt in enumerate(points):
        
        mirrored = np.array(([pt[0], pt[1], -pt[2]]))
        if mirrored not in points:
            print("The following vertex is not mirrored: \n")
            print(pt)

    return

def verify_xy_symmetry(points):

    # Iterate one point at a time
    for i, pt in enumerate(points):
        check = False
        
        mirrored = np.array(([pt[0], -pt[1], pt[2]]))
        if mirrored not in points:
            print("The following vertex is not mirrored: \n")
            print(pt)

    return



# Main
# Read in VTK file
filename = 'meshes/tapered_wing_test.vtk'
mesh = pv.read(filename)
vertices = mesh.points

# Check symmetry of xz planes
verify_xz_symmetry(vertices)

# Check symmetry of xy planes
verify_xy_symmetry(vertices)

print('Mirroring verification complete')