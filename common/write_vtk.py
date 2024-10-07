import vtk
import numpy as np


def add_data_to_vtk(existing_vtk_file, processed_CFx, processed_CFy, processed_CFz, output_vtk_file):
    # Read the existing VTK file
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(existing_vtk_file)
    reader.Update()

    polydata = reader.GetOutput()

    # Create new vtkFloatArrays to hold the vector data
    d_CFx_array = vtk.vtkFloatArray()
    d_CFx_array.SetNumberOfComponents(3)
    d_CFx_array.SetName("CFx Vectors")

    d_CFy_array = vtk.vtkFloatArray()
    d_CFy_array.SetNumberOfComponents(3)
    d_CFy_array.SetName("CFy Vectors")

    d_CFz_array = vtk.vtkFloatArray()
    d_CFz_array.SetNumberOfComponents(3)
    d_CFz_array.SetName("CFz Vectors")

    # Insert tuples for each vector component
    for vec in processed_CFx:
        d_CFx_array.InsertNextTuple3(*vec)  # Unpack tuple directly into the method

    for vec in processed_CFy:
        d_CFy_array.InsertNextTuple3(*vec)

    for vec in processed_CFz:
        d_CFz_array.InsertNextTuple3(*vec)

    # Add the vector arrays to the polydata
    polydata.GetCellData().AddArray(d_CFx_array)
    polydata.GetCellData().AddArray(d_CFy_array)
    polydata.GetCellData().AddArray(d_CFz_array)

    # Write the updated polydata to a new VTK file
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(output_vtk_file)
    writer.SetInputData(polydata)
    writer.Write()
    print(f"Updated VTK file written to {output_vtk_file}")











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