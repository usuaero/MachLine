import paraview.simple as pvs
import numpy as np

AoA_folder = ['0', '2.1', '4.2', '4.3', '6.2', '6.3', '6.4', '8.3', '8.4', '8.5', '8.6', '10.4', '10.5', '10.6', '10.7']
locations = {
    '4.1': ['0', '2.1', '4.2','6.2','8.3','10.4'], #angles of attack for 4.1% semispan
    '8.2': ['0', '2.1', '4.2','6.2','8.3','10.4'], #angles of attack for 8.2% semispan
    '16.3': ['0', '2.1', '4.2','6.2','8.3','10.4'], #angles of attack for 16.3% semispan
    '24.5': ['0', '2.1', '4.2','6.2','8.3','10.4'], #angles of attack for 24.5% semispan
    '36.7': ['0', '2.1', '4.2','6.3','8.4','10.5'], #angles of attack for 36.7% semispan
    '51.0': ['0', '2.1', '4.2','6.3','8.4','10.5'], #angles of attack for 51.0% semispan
    '65.3': ['0', '2.1', '4.2','6.3','8.5','10.6'], #angles of attack for 65.3% semispan
    '89.8': ['0', '4.3', '10.7'], #angles of attack for 89.8% semispan
    '94.9': ['0', '2.1', '4.3','6.4','8.6','10.7'] #angles of attack for 94.9% semispan
}

semispan_xy_loc = [2, 4, 7.99, 12, 17.98, 24.99, 32, 44, 46.5]
points_to_process = ['C_p', 'sigma', 'v']

Nodes = ['10', '17', '25', '45']

#Iterate over Node Density Files
for Node in Nodes:

    # Iterate over all angles of attack
    for AoA in AoA_folder:

        # Identify the filename and filepath of the vtk file being processed
        vtk_file_location = 'MachLine_Results/' + AoA + '_degrees_AoA/half_wing_A_' + Node + '_nodes_' + AoA + '_deg_AoA.vtk'
        PassedName= 'half_wing_A_' + Node + '_nodes_' + AoA + '_deg_AoA.vtk'

        # Load vtk file into paraview
        data_reader = pvs.LegacyVTKReader(registrationName=PassedName, FileNames= vtk_file_location)

        # create a new 'Cell Data to Point Data'
        cellDatatoPointData1 = pvs.CellDatatoPointData(registrationName='CellDatatoPointData1', Input=data_reader)
        cellDatatoPointData1.CellDataArraytoprocess = points_to_process

        # create a new 'Plot On Intersection Curves'
        plotOnIntersectionCurves1 = pvs.PlotOnIntersectionCurves(registrationName='PlotOnIntersectionCurves1', Input=cellDatatoPointData1)
        plotOnIntersectionCurves1.SliceType = 'Plane'



        # Iterate over percent semispan locations     
        for i, percent_semispan in enumerate(locations):

            # Verify that current AoA is applicable at selected semispan location
            if AoA in locations[percent_semispan]:

                # Properties modified on plotter.SliceType
                set_origin = [semispan_xy_loc[i], semispan_xy_loc[i], 0]
                set_normal = [0,1,0]
                plotOnIntersectionCurves1.SliceType.Origin = set_origin
                plotOnIntersectionCurves1.SliceType.Normal = set_normal

                # Create a new 'Line Chart View'
                lineChartView1 = pvs.CreateView('XYChartView')

                # show data in view
                plotOnIntersectionCurves1Display = pvs.Show(plotOnIntersectionCurves1, lineChartView1, 'XYChartRepresentation')

                # update the view to ensure updated data information
                lineChartView1.Update()

                # Properties modified on plotOnIntersectionCurves1Display
                plotOnIntersectionCurves1Display.XArrayName = 'Points_X'

                # update the view to ensure updated data information
                lineChartView1.Update()

                # save data
                save_location = percent_semispan + '_percent_semispan/' + AoA + '_degrees_AoA/' + percent_semispan + '_percent_semispan_' + Node + '_nodes_results.csv'
                pvs.SaveData(save_location, proxy=plotOnIntersectionCurves1, PointDataArrays=points_to_process, FieldAssociation='Point Data', Precision=12)

                

                