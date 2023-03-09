import paraview.simple as pvs
import numpy as np
import json
import os

# Import geometry and processing locations from the input file
filename = "Swept_half_wing_conditions_input.json"

json_string = open(filename).read()
json_vals = json.loads(json_string)

AoA_folder = json_vals["geometry"]["AoA list"]
locations = json_vals["geometry"]["semispan locations and AoA"]

semispan_xy_loc = json_vals["geometry"]["semispan_xy_loc"]
points_to_process = json_vals["plots"]["points to process"]

Nodes = json_vals["geometry"]["nodes"]

# Import formulation 
formulation_inputs = json_vals["solver"]["formulation"]

print("\nData processing in progress...\n")

# Iterate over formulations
for form in formulation_inputs:

    if form == "source-free":
        formulation_adjusted = "source_free"
    else:
        formulation_adjusted = form

    #Iterate over Node Density Files
    for Node in Nodes:

        # Iterate over all angles of attack
        for AoA in AoA_folder:

            # Identify the filename and filepath of the vtk file being processed
            PassedName= "half_wing_A_" + Node + "_nodes_" + AoA + "_deg_AoA_" + formulation_adjusted + "_formulation.vtk"
            vtk_file_location = 'MachLine_Results/' + AoA + '_degrees_AoA/' + PassedName

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
                    set_origin = [semispan_xy_loc[percent_semispan], semispan_xy_loc[percent_semispan], 0]
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
                    save_folder = percent_semispan + '_percent_semispan/' + AoA + '_degrees_AoA/'
                    if not os.path.exists(save_folder):
                        os.makedirs(save_folder)
                    save_location = save_folder + percent_semispan + '_percent_semispan_' + Node + '_nodes_results_' + AoA + '_deg_AoA_' + formulation_adjusted + '_formulation.csv'
                    pvs.SaveData(save_location, proxy=plotOnIntersectionCurves1, PointDataArrays=points_to_process, FieldAssociation='Point Data', Precision=12)

                    

print("Data processing complete.\n")