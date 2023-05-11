import sys
import copy
import numpy as np

if __name__=="__main__":

    # Get file
    mesh_file = sys.argv[-1]

    # Open
    with open(mesh_file, 'r') as mesh_file_handle:

        # get lines
        lines = mesh_file_handle.readlines()

    # Find where panels start
    start = 0
    while "POLYGONS" not in lines[start]:
        start += 1
    start += 1

    # Find where panels end
    end = start
    while "CELL_DATA" not in lines[end]:
        end += 1

    # Get panel indices
    panel_data = np.zeros((3,end-start), dtype=int)
    for i in range(end-start):
        panel_data[:,i] = np.array(lines[i+start].split()[1:], dtype=int)

    # Rewrite
    with open(mesh_file.replace(".vtk", "_fixed.vtk"), 'w') as mesh_file_handle:

        # Write out original starting lines
        for i in range(start):
            print(lines[i], file=mesh_file_handle, end='')

        # Write out fixed panel data
        for i in range(end-start):
            line = "3 " + " ".join([str(x) for x in panel_data[::-1,i]])
            print(line, file=mesh_file_handle)