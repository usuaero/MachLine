import sys
import vtk

print(vtk)

def convert_stl_to_vtk(stl_file):
    a = vtk.vtkSTLReader()
    a.SetFileName(stl_file)
    a.Update()
    a = a.GetOutput()

    new_file = stl_file.replace('.stl', '.vtk')
    b = vtk.vtkPolyDataWriter()
    b.SetFileName(new_file)
    b.SetInputData(a)
    b.Update()

if __name__=="__main__":

    for mesh_file in sys.argv:
        if ".stl" not in mesh_file:
            continue
        convert_stl_to_vtk(mesh_file)