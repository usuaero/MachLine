import vtk
import numpy as np

def read_vtk_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    return lines

def write_vtk_file(filename, lines):
    with open(filename, 'w') as file:
        file.writelines(lines)

def add_vector_data_to_vtk(lines, vector_data_sets):
    # Update points (assuming new points are to be added)
    num_points = len(vector_data_sets[0])  # Assuming all sets have the same number of points

    updated_lines = lines

    updated_lines.append(f'POINT_DATA                  {num_points}\n')
    updated_lines.append(f'VECTORS d_CFx_central_diff float\n')
    for vector in vector_data_sets[0]:
            updated_lines.append(f'{vector[0]: .12E} {vector[1]: .12E} {vector[2]: .12E}\n')

    updated_lines.append(f'VECTORS d_CFy_central_diff float\n')
    for vector in vector_data_sets[1]:
            updated_lines.append(f'{vector[0]: .12E} {vector[1]: .12E} {vector[2]: .12E}\n')

    updated_lines.append(f'VECTORS d_CFz_central_diff float\n')
    for vector in vector_data_sets[2]:
            updated_lines.append(f'{vector[0]: .12E} {vector[1]: .12E} {vector[2]: .12E}\n')

    # # Add vector datasets
    # for i, vector_set in enumerate(vector_data_sets):
    #     updated_lines.append(f'VECTORS vector_set_{i+1} float\n')  # New data header for each vector set
    #     for vector in vector_set:
    #         updated_lines.append(f'{vector[0]: .12E} {vector[1]: .12E} {vector[2]: .12E}\n')

    return updated_lines


# class VTK_Poly_Data_Writer:

#     # initalize with the desired vtk file name and location and the set of geometric points
#     def __init__(self, filename, points):

#         # store the file name
#         self.filename = filename

#         # pull out the point data
#         for i in range(0, len(points), 3):
#             self.points = []
        
#         # initialize a list of vector sets
#         self.vector_sets = {}


#     # add a vector the to data that will be written in a vtk
#     def add_vector_data(self, vector_name, vector_data):

#         # pull out each vector that will be associated with a point
#         for i in range(0, len(vector_data), 3):
#             vtk_vectors = []

#         # add this vector set to the vector set list attribute
#         self.vector_sets[vector_name] = vtk_vectors


#     # write the vtk file
#     def write_vtk(self):

#         # initialize VTK data structure for data associated with points
#         point_data = pyvtk.PointData()

#         # add each vector set and its name to point_data
#         for name, vector_field in self.vector_sets.items():
#             point_data.append(pyvtk.Vectors(vector_field, name = name))

#         # create the polydata type or something like that
#         vtk_data = pyvtk.VtkData(pyvtk.PolyData(points = self.points), point_data)


#         # write to file
#         vtk_data.tofile(self.filename)