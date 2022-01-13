import paraview.simple as pvs

AoA_folder = ['0', '2_1', '4_2', '4_3', '6_2', '6_3', '6_4', '8_3', '8_4', '8_5', '8_6', '10_4', '10_5', '10_6', '10_7']

Nodes = ['10', '17', '25', '45']

for i in range(len(AoA_folder)):
    for j in range(len(Nodes)):
        vtk_file_location = 'MachLine_Results/' + AoA_folder[i] + '_degrees_AoA/half_wing_A_'
        if AoA_folder[i] == '0':
            vtk_file_location += 'updated_mesh_' + Nodes[j] + '_nodes.vtk'

        else:
            vtk_file_location += Nodes[j] + '_nodes_' + AoA_folder[i] + '_deg_AoA.vtk'

        # Load vtk file into paraview
        data_reader = pvs.OpenDataFile(vtk_file_location)