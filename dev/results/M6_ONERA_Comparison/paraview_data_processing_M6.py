import paraview.simple as pvs
import numpy as np
import json
from sys import exit

# Import geometry and processing locations from the input file
filename = "M6_Onera_input_settings.json"

json_string = open(filename).read()
json_vals = json.loads(json_string)

AoA_folder = json_vals["geometry"]["AoA list"]
mach_numbers = json_vals["flow conditions"]["mach_numbers"]
length = json_vals["geometry"]["semispan_length"]
sweep_angle = json_vals["geometry"]["sweep_angle[deg]"] * np.pi / 180
locations = json_vals["geometry"]["semispan locations and AoA"]
semispan_y_loc = json_vals["geometry"]["semispan_y_loc"]
datasets = json_vals["geometry"]["datasets"]

# Import solver settings
formulation_inputs = json_vals["solver"]["formulation"]
calculations = json_vals["solver"]["calculation_type"]

print("\nData processing in progress...\n")

# Iterate over formulations
for form in formulation_inputs:

    if form == "source-free":
        formulation_adjusted = "source_free"
    else:
        formulation_adjusted = form

    # Iterate over data file types
    for data in datasets:

        # Determine points to process based on what resource the data was pulled from
        if data == "FlightStream":
            points_to_process = ["Cp (1)"]
        else:
            points_to_process = json_vals["plots"]["points_to_process"]

        # Iterate over calculation types, direct subsonic calculations, or post-processing pressure corrections
        for calc in calculations:

            # Iterate over all angles of attack
            for k, AoA in enumerate(AoA_folder):

                # Identify the filename and filepath of the vtk file being processed
                if data == "MachLine":
                    PassedName= f"{data}_mach_" + mach_numbers[k] + "_AoA_" + AoA + "_" + formulation_adjusted + ".vtk"
                    vtk_file_location = f'M6_meshes/{data}/{calc}/' + PassedName
                elif data == "FlightStream":
                    PassedName= f"{data}_mach_" + mach_numbers[k] + "_AoA_" + AoA + ".vtk"
                    vtk_file_location = f'M6_meshes/{data}/' + PassedName
                else:
                    print("Invalid Data Type Selected. Quitting...")
                    exit()


                # Load vtk file into paraview
                data_reader = pvs.LegacyVTKReader(registrationName=PassedName, FileNames= vtk_file_location)

                # create a new 'Cell Data to Point Data' unless FlightStream is used
                if data != "FlightStream":
                    cellDatatoPointData1 = pvs.CellDatatoPointData(registrationName='CellDatatoPointData1', Input=data_reader)
                    cellDatatoPointData1.CellDataArraytoprocess = points_to_process
                    # create a new 'Plot On Intersection Curves'
                    plotOnIntersectionCurves1 = pvs.PlotOnIntersectionCurves(registrationName='PlotOnIntersectionCurves1', Input=cellDatatoPointData1)
                else:
                    passed_data = pvs.FindSource(PassedName)
                    plotOnIntersectionCurves1 = pvs.PlotOnIntersectionCurves(registrationName='PlotOnIntersectionCurves1', Input=passed_data)

                plotOnIntersectionCurves1.SliceType = 'Plane'

                # Iterate over percent semispan locations     
                for i, percent_semispan in enumerate(locations):

                    # Verify that current AoA is applicable at selected semispan location
                    if AoA in locations[percent_semispan]:

                        # Calculate x location of LE based on sweep angle
                        x_location = float(percent_semispan) * np.sin(sweep_angle)

                        # Properties modified on plotter.SliceType
                        set_origin = [x_location, semispan_y_loc[percent_semispan], 0]
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
                        if data == "FlightStream":
                            plotOnIntersectionCurves1Display.SeriesVisibility = points_to_process
                        # update the view to ensure updated data information
                        lineChartView1.Update()

                        # save data
                        if data == "MachLine":
                            save_location = f"{data}_Data/{calc}/M6_{data}_data_{percent_semispan}_percent_semispan_mach_{mach_numbers[k]}_AoA_{AoA}"'.csv'
                        else:
                            save_location = f"{data}_Data/M6_{data}_data_{percent_semispan}_percent_semispan_mach_{mach_numbers[k]}_AoA_{AoA}"'.csv'
                        pvs.SaveData(save_location, proxy=plotOnIntersectionCurves1, PointDataArrays=points_to_process, FieldAssociation='Point Data', Precision=12)
                
                    

print("Data processing complete.\n")