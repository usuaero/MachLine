#This file is to take .csv or .txt files of lift distributions and plot
#them allowing comparison to experimental results

import os
import json
import numpy as np
import matplotlib.pyplot as plt

class Swept_Plotting:

    def __init__(self, Geometry_Data, data_type, y_locations, process_points_dict, calculation_type):

        self.aoa_for_locations = Geometry_Data
        self.semispan_y_loc = y_locations
        self.data_type = data_type
        self.points_dict = process_points_dict
        self.pressure_rule = calculation_type


    def plot_machline_pressures(self, file_name, pressure_rule):

        # Create a lists to distinguish direct subsonic pressure results and pressure correction resutls
        pressure_corrections = ["C_p_inc", "C_p_PG", "C_p_KT", "C_p_L"]
        direct_subsonic = ["C_p_ise", "C_p_2nd", "C_p_lin", "C_p_sln"]

        # Create an updated dictionary with correct pressure coefficients based on rule type
        dict_copy = self.points_dict.copy()
        
        # Iterate over keys in dictionary
        for key in self.points_dict:
            if pressure_rule == "pressure_corrections":
                if key in direct_subsonic:
                    del dict_copy[key]
            if pressure_rule == "subsonic_calcs":
                if key in pressure_corrections:
                    del dict_copy[key]
 
        # Get data in string format
        data_str = np.genfromtxt(file_name, delimiter=",", skip_header=0, dtype=str)

        # Remove multiple quotes
        data_str = np.char.replace(data_str,'"',"")

        # Get data in floating-point format
        data_flp = np.genfromtxt(file_name, delimiter=",", skip_header=1, dtype=float)

        # Iterate over points to process
        for j, point in enumerate(dict_copy):

            # Verify requested point is in the dataset
            if point not in data_str[0,:]:
                print('!!! {0} is not found in {1}. Quitting...'.format(point, file_name))
                exit()

            # Iterate over file columns searching for the x and Cp columns
            for i, val in enumerate(data_str[0,:]):
                if val == 'Points:0':
                    x_loc = i
                if val == point:
                    point_loc = i
                
            point_data = data_flp[:,point_loc]
            x = data_flp[:,x_loc]
            
            # Calculate local chord length and LE location for each section analyzed
            Chord_Length = max(x) - min(x)
            LE_loc = min(x)
            x_axis = (x-LE_loc)/Chord_Length
            
            # Black and white lines with line types and varying opacity
            if pressure_rule == "subsonic_calcs":
                plt.plot(x_axis, point_data, 'k-', label="Direct")
            elif pressure_rule == "pressure_corrections":
                plt.plot(x_axis, point_data, 'k--', label="P-G Corr.")


    def plot(self, mach_number):

        # Pull In Results Depending on Location Input
        # In addition to defining leading edge location and respective labels for plot
        for percent_semispan in self.aoa_for_locations.keys():

            # Define section number based on geometry of wing
            if percent_semispan == "20":
                section = 0
            elif percent_semispan == "44":
                section = 1
            elif percent_semispan == "65":
                section = 2
            elif percent_semispan == "80":
                section = 3
            elif percent_semispan == "90":
                section = 4
            elif percent_semispan == "95":
                section = 5
            elif percent_semispan == "99":
                section = 6
            else:
                print("Invalid percent semispan. Quitting...")
                exit()

            # Establish clean list for saving and combining produced plots
            list_of_files =[]
                
            # Iterate over angles of attack for each semispan
            for j, AoA in enumerate(self.aoa_for_locations[percent_semispan]):

                # Pull in experimental results
                file_experimental = "studies/M6_ONERA_Comparison/Experimental_Data/M6_Data_mach_{0}_AoA_{1}.csv".format(mach_number[j], AoA)
                exp_data = np.genfromtxt(file_experimental, delimiter=",", skip_header=2, dtype=float, usecols=(section*2, section*2+1))
                plt.figure()
                plt.plot(exp_data[:,0], exp_data[:, 1], "ko", label="Experiment", markersize=3)

                # Iterate over data sources
                for data_type in self.data_type:
                    if data_type == "MachLine":

                        # Iterate over pressure rule
                        for rule in self.pressure_rule:
                            file = "studies/M6_ONERA_Comparison/csv_file_data/{0}/M6_{1}_data_{2}_percent_semispan_mach_{3}_AoA_{4}.csv".format(rule, data_type, percent_semispan, mach_number[j], AoA)
                            
                            # Call plotting function
                            self.plot_machline_pressures(file, rule)
                    else:
                        file = f"studies/M6_ONERA_Comparison/{data_type}_Data/M6_{data_type}_data_{percent_semispan}_percent_semispan_mach_{mach_number[j]}_AoA_{AoA}.csv"
                    
                        # Call plotting function
                        self.plot_machline_pressures(file, data_type)
                    
                # Invert the Y-axis
                plt.gca().invert_yaxis()

                # Identify labels and notes to include on plots
                title = f"Pressure Distribution at section {section + 1}, {percent_semispan} percent semispan"

                # Add labels
                plt.xlabel("$x/c$")
                plt.ylabel("$C_P$")
                plt.legend()
                
                # Only include title and notes if the plots aren't going to be saved
                if (not json_vals["plots"]["save plots"]) or json_vals["plots"]["combine pdfs"]:
               
                    #Identify angle of atack and mach number notes to add to plot title
                    y = '$\\alpha$='
                    AoA_Notes = y + AoA + '$^{\circ}$'
                    M_notes = '$M_{\infty}$ = ' + mach_number[j] + ","

                    plt.title(title)
                    plt.figtext(0.675, 0.01, M_notes)
                    plt.figtext(0.825, 0.015, AoA_Notes)
                
                # Save the figure in appropriate location
                if json_vals["plots"]["save plots"]:

                    # Initialize plot directory
                    plot_dir = "studies/M6_ONERA_Comparison/Plots/"
                    if not os.path.exists(plot_dir):
                        os.makedirs(plot_dir)

                    # Save file
                    filename = plot_dir + "M6_section_{0}_mach_{1}_aoa_{2}.".format(section+1, mach_number[j], AoA) + json_vals["plots"]["save plot type"]
                    list_of_files.append(filename)
                    plt.savefig(filename)
                    plt.close()

                else:
                    plt.show()


if __name__=="__main__":

    # Main script
    input_file = "studies/M6_ONERA_Comparison/M6_Onera_input_settings.json"
    json_string = open(input_file).read()
    json_vals = json.loads(json_string)

    # Identify geometry and flow conditions
    aoa_for_locations = json_vals["geometry"]["semispan locations and AoA"]
    semispan_y_loc = json_vals["geometry"]["semispan_y_loc"]
    datasets = json_vals["geometry"]["datasets"]
    mach_numbers = json_vals["flow conditions"]["mach_numbers"]

    # Identify solver and plotting information
    formulation = json_vals["solver"]["formulation"]
    subsonic_calc_type = json_vals["solver"]["calculation_type"]
    points_to_process = json_vals["plots"]["points_to_process"]

    # Format the points to process for plotting purposes
    temp = np.copy(points_to_process)
    temp = np.char.replace(temp, "_", "",1)
    new_temp = []
    for name in temp:
        for j, letter in enumerate(name):
            if letter == "_":
                spot = j
        new_name = name[:spot+1] + "{" + name[spot+1:] + "}"

        # Account for changing isentropic to direct (ise -> dir)
        if 'ise' in new_name:
            new_name = np.char.replace(new_name, "ise", "dir",1)
        new_temp.append(new_name)

    # Create a dictionary to contain the legend names properly formatted for the requested parameters
    labels_dict = dict(zip(points_to_process, new_temp))

    # Iterate over solver formulations from input file

    # Proccesses either one or all semispan aoa_for_locations based on input file
    if json_vals["plots"]["process_all"]:

        Swept_Plotting(aoa_for_locations, datasets, semispan_y_loc, labels_dict, subsonic_calc_type).plot(mach_numbers)

    else:
        x = input("Enter percent semispan for analysis results. Options are 20, 44, 65, 80, 90, 95, 99:   ",)
        if x in aoa_for_locations.keys():
            Specific_Semispan = {x: aoa_for_locations[x]}

            Swept_Plotting(Specific_Semispan, datasets, semispan_y_loc, labels_dict, subsonic_calc_type).plot(mach_numbers)

        else:
            print("\n****************\nInvalid Entry. Please run script again.\n****************\n")
            quit()

    # Print exit statement to verify completion of process
    print()
    print("Pressure plotting executed successfully. \n")