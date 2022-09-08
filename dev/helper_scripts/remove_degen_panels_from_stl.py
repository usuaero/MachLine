import sys
import numpy as np

if __name__=="__main__":

    # Get mesh file
    mesh_file = sys.argv[-1]

    # Check for an actual file
    if mesh_file == "remove_degen_panels_from_stl.py":
        raise IOError

    # Read in lines
    with open(mesh_file, 'r') as mesh_file_handle:
        lines = mesh_file_handle.readlines()

    # Number of panels
    N_panels = (len(lines)-2)//7
    is_degen = np.zeros(N_panels, dtype=bool)

    # Loop through panels
    for i in range(N_panels):

        # Get vertices
        vert1 = np.array([float(x) for x in lines[i*7+3].split()[1:]])
        vert2 = np.array([float(x) for x in lines[i*7+4].split()[1:]])
        vert3 = np.array([float(x) for x in lines[i*7+5].split()[1:]])

        # Area
        A = 0.5*np.linalg.norm(np.cross(vert2-vert1, vert3-vert2))

        # Check
        if A < 1e-12:
            print("Panel {0} is degenerate.".format(i))
            is_degen[i] = True
        else:
            is_degen[i] = False

    # Write back out non-degenerate panels
    with open(mesh_file.replace(".stl", "_fixed.stl"), 'w') as mesh_file_handle:

        # Write header
        mesh_file_handle.write(lines[0])

        # Write non-degenerate panels
        for i in range(N_panels):
            if not is_degen[i]:
                mesh_file_handle.writelines(lines[i*7+1:i*7+8])

        # Write footer
        mesh_file_handle.write(lines[-1])