import numpy as np
import vtk

def read_vtk(vtk_filename):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(vtk_filename)
    reader.Update()

    data = reader.GetOutput()
    print(data)
    cfx = data.GetPointData().GetVectors('CF_x_sensitivity_vectors')
    cfy = data.GetPointData().GetVectors('CF_y_sensitivity_vectors')
    cfz = data.GetPointData().GetVectors('CF_z_sensitivity_vectors')
    normals = data.GetPointData().GetVectors('vertex_outward_normal_vectors')
    # freestream = data.GetPointData().GetScalars('freestream')
    
    # print("Freestream vector:", freestream)
    # if freestream:
    #     freestream_vector = freestream.GetTuple(0)
    # else:
    #     print("freestream vector not found in vtk")
    freestream_vector = [100,0,10.0,10.0]
    num_points = data.GetNumberOfPoints()
    cfx_vectors = np.array([cfx.GetTuple3(i) for i in range(num_points)])
    cfy_vectors = np.array([cfy.GetTuple3(i) for i in range(num_points)])
    cfz_vectors = np.array([cfz.GetTuple3(i) for i in range(num_points)])
    point_normals = np.array([normals.GetTuple3(i) for i in range(num_points)])
    # freestream_vector = np.array(freestream.GetTuple3(0))
    
    return num_points, cfx_vectors, cfy_vectors, cfz_vectors, point_normals, freestream_vector, reader


def update_vtk(vtk_filename, CL_adjoint, CD_adjoint, CS_adjoint, reader):
    
    data = reader.GetOutput()

    # write in CL sensitivity data
    CL_data = vtk.vtkDoubleArray()
    CL_data.SetNumberOfComponents(1)
    CL_data.SetName("CL_norm_sensitivities")
    for value in CL_adjoint:
        CL_data.InsertNextValue(value)
    data.GetPointData().AddArray(CL_data)

    # write in CD sensitivity data
    CD_data = vtk.vtkDoubleArray()
    CD_data.SetNumberOfComponents(1)
    CD_data.SetName("CD_norm_sensitivities")
    for value in CD_adjoint:
        CD_data.InsertNextValue(value)
    data.GetPointData().AddArray(CD_data)

    # write in CS sensitivity data
    CS_data = vtk.vtkDoubleArray()
    CS_data.SetNumberOfComponents(1)
    CS_data.SetName("CS_norm_sensitivities")
    for value in CS_adjoint:
        CS_data.InsertNextValue(value)
    data.GetPointData().AddArray(CS_data)

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(vtk_filename)
    writer.SetInputData(data)
    writer.Write()
    
if __name__ == "__main__":
    vtk_filename = 'dev/results/adjoint/octa_mesh_body_results.vtk'

    num_points, cfx_vectors, cfy_vectors, cfz_vectors, point_normals, freestream_vector, reader = read_vtk(vtk_filename)
    u = freestream_vector[0]
    v = freestream_vector[1]
    w = freestream_vector[2]

    alpha = np.arctan2(u,w)
    beta = np.arcsin(v/np.sqrt(u*u + v*v + w*w))
    
    CL_sensitivities = cfz_vectors*np.cos(alpha) - cfx_vectors*np.sin(alpha) 
    CD_sensitivities = cfx_vectors*np.cos(alpha)*np.cos(beta) - cfy_vectors*np.sin(beta) + cfz_vectors*np.sin(alpha)*np.cos(beta) 
    CS_sensitivities = cfx_vectors*np.cos(alpha)*np.sin(beta) + cfy_vectors*np.cos(beta) + cfz_vectors*np.sin(alpha)*np.sin(beta)

    CL_norm_sensitivitites=np.zeros(num_points) 
    CD_norm_sensitivitites=np.zeros(num_points) 
    CS_norm_sensitivitites=np.zeros(num_points)
    
    for i in range(num_points):
        CL_norm_sensitivitites[i] = np.dot(CL_sensitivities[i,:], point_normals[i,:])
        CD_norm_sensitivitites[i] = np.dot(CD_sensitivities[i,:], point_normals[i,:])
        CS_norm_sensitivitites[i] = np.dot(CS_sensitivities[i,:], point_normals[i,:])

    update_vtk(vtk_filename, CL_norm_sensitivitites, CD_norm_sensitivitites, CS_norm_sensitivitites, reader)