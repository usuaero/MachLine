import copy
import sys
import pypan as pp
import numpy as np

if __name__=="__main__":

    for mesh_file in sys.argv:
        if ".stl" not in mesh_file:
            continue

        my_mesh = pp.Mesh(name="", mesh_file=mesh_file, mesh_file_type="STL", kutta_angle=90.0, verbose=True)
        my_mesh.export_vtk(mesh_file.replace(".stl", ".vtk"))
