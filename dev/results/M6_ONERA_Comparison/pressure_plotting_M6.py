#This file is to take .csv or .txt files of lift distributions and plot
#them allowing comparison to experimental results

import numpy as np
import matplotlib.pyplot as plt
import json
from PyPDF2 import PdfFileMerger
from sklearn.metrics import precision_recall_curve

class Swept_Plotting:

    def __init__(self, Geometry_Data, data_type, y_locations, process_points_dict, calculation_type):

        self.locations = Geometry_Data
        self.semispan_y_loc = y_locations
        self.data_type = data_type
        self.points_dict = process_points_dict
        self.calc_type = calculation_type

    def Pressure_Plot(self, file_name, calc_type):

        # Create a lists to distinguish direct subsonic pressure results and pressure correction resutls
        pressure_corrections = ["C_p_inc", "C_p_PG", "C_p_KT", "C_p_L"]
        direct_subsonic = ["C_p_ise", "C_p_2nd", "C_p_lin", "C_p_sln"]

        # Create an updated dictionary with correct pressure coefficients based on calc type
        dict_copy = self.points_dict.copy()
        
        # Iterate over keys in dictionary
        for key in self.points_dict:
            if calc_type == "pressure_corrections":
                if key in direct_subsonic:
                    del dict_copy[key]
            if calc_type == "subsonic_calcs":
                if key in pressure_corrections:
                    del dict_copy[key]
 
        Data_complete = np.genfromtxt(file_name, delimiter=",", skip_header=0, dtype=str)
        # Reformate Data_complete to remove multiple quotes
        Data_complete = np.char.replace(Data_complete,'"',"")

        Data = np.genfromtxt(file_name, delimiter=",", skip_header=1, dtype=float)

        # Create lists to pull from for plotting differentiation
        if calc_type == "pressure_corrections":
            opacity = [1, 0.5, 0.3,]
            shape = ["v", "^", "<"]
            linetype = ["--"]
        elif calc_type == "subsonic_calcs":
            opacity = [1, 0.5, 0.3,]
            shape = ["o", "s", "p"]
            linetype= ["-"]
            # colors = ["g", "b", "r", "m", "y"]

        # Iterate over points to process
        for j, point in enumerate(dict_copy):
            # Verify requested point is in the dataset
            if point not in Data_complete[0,:]:
                print(f'***{point} is not found in {file_name}. Quitting')
                exit()

            # Iterate over file columns searching for the x and Cp columns
            for i, val in enumerate(Data_complete[0,:]):
                if val == 'Points:0':
                    x_loc = i
                if val == point:
                    point_loc = i
                
            point_data = Data[:,point_loc]
            x = Data[:,x_loc]
            
            # Calculate local chord length and LE location for each section analyzed
            Chord_Length = max(x) - min(x)
            LE_loc = min(x)
            x_axis = (x-LE_loc)/Chord_Length
            
            # Black and white lines with line types and varying opacity
            plt.plot(x_axis, point_data, label=f"${dict_copy[point]}$", linestyle=linetype[j], alpha=opacity[j], color="k")
            
            # Various colors
            # plt.plot((x-LE_xy_loc)/Chord_Length,point_data, "o", label=f"${dict_copy[point]}$", color=colors[j])
            
            # Markers alone without lines for all data
            # plt.plot(x_axis,point_data, label=f"${dict_copy[point]}$", marker=shape[j], alpha=opacity[j], color="k", linestyle="none")

        return



    def get_data(self, mach_number):


        # Pull In Results Depending on Location Input
        # In addition to defining leading edge location and respective labels for plot
        for percent_semispan in self.locations:

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
            for j,AoA in enumerate(self.locations[percent_semispan]):
                # Pull in experimental results
                file_experimental = f"Experimental_Data/M6_Data_mach_{mach_number[j]}_AoA_{AoA}.csv"
            
                Experimental = np.genfromtxt(file_experimental, delimiter=",", skip_header=2, dtype=float, usecols=(section*2, section*2+1))
                # Iterate over data sources
                for data_type in self.data_type:
                    if data_type == "MachLine":

                        # Iterate over calculation type
                        for calc in self.calc_type:
                            file = f"{data_type}_Data/{calc}/M6_{data_type}_data_{percent_semispan}_percent_semispan_mach_{mach_number[j]}_AoA_{AoA}.csv"
                            
                            # Call plotting function
                            self.Pressure_Plot(file, calc)
                    else:
                        file = f"{data_type}_Data/M6_{data_type}_data_{percent_semispan}_percent_semispan_mach_{mach_number[j]}_AoA_{AoA}.csv"
                    
                        # Call plotting function
                        self.Pressure_Plot(file, data_type)
                    
                # Invert the Y-axis
                plt.gca().invert_yaxis()


                #=== Plot Details ===
                #Identify labels and notes to include on plots
                title = f"Pressure Distribution at section {section + 1}, {percent_semispan} percent semispan"
               
                #Identify angle of atack and mach number notes to add to plot title
                y = r'$\alpha$='
                AoA_Notes = y + AoA + r'$^{\circ}$'
                M_notes = r'$M_{\infty}$ = ' + mach_number[j] + ","

                #Pull in experimental results for comparison
                plt.plot(Experimental[:,0], Experimental[:, 1], "o", label="Experimental", color="k", fillstyle="full")
                # plt.plot(Experimental[:,0], Experimental[:, 1], label="Experimental", color="k")

                #Plot the figure containing all curves
                xlabel = 'x/c'
                plt.xlabel(xlabel)
                plt.ylabel(r"$C_p$")
                plt.legend()
                
                # Only include title and notes if the plots aren't going to be saved
                if (not json_vals["plots"]["save plots"]) or json_vals["plots"]["combine pdfs"]:
                    plt.title(title)
                    plt.figtext(0.675, 0.01, M_notes)
                    plt.figtext(0.825, 0.015, AoA_Notes)
                # complete_title = title + ", " + AoA_Notes
                # Sub_Note = "*" + formulation + " form."
                # Identify location of formulation type on plot
                # plt.title(complete_title)
                
                #Save the figure in appropriate location and with the correct size
                if json_vals["plots"]["save plots"]:
                    filename = f"Plots/M6_Onera_comparison_mach_{mach_number[j]}_section_{section+1}." + json_vals["plots"]["save plot type"]
                    list_of_files.append(filename)
                    plt.savefig(filename)
                plt.show()

            # Combine plots into one pdf for easy review
            if json_vals["plots"]["save plot type"] == "pdf" and json_vals["plots"]["combine pdfs"]:

                # Combine files in Combined Plots Folder
                merger = PdfFileMerger()
                Combined_name = f"Plots/Combined_Plots/M6_Onera_comparison_{percent_semispan}_percent_semispan." + json_vals["plots"]["save plot type"]

                for name in list_of_files:
                    merger.append(name)

                merger.write(Combined_name)
                merger.close()

# Main script
inputfile = "M6_Onera_input_settings.json"
json_string = open(inputfile).read()
json_vals = json.loads(json_string)

# Identify geometry and flow conditions
locations = json_vals["geometry"]["semispan locations and AoA"]
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
    new_temp.append(new_name)

# Create a dictionary to contain the legend names properly formatted for the requested parameters
labels_dict = dict(zip(points_to_process, new_temp))

# Iterate over solver formulations from input file

    # Proccesses either one or all semispan locations based on input file
if json_vals["plots"]["process_all"]:

    Swept_Plotting(locations, datasets, semispan_y_loc, labels_dict, subsonic_calc_type).get_data(mach_numbers)

else:
    x = input("Enter percent semispan for analysis results. Options are 20, 44, 65, 80, 90, 95, 99:   ",)
    if x in locations.keys():
        Specific_Semispan = {x: locations[x]}

        Swept_Plotting(Specific_Semispan, datasets, semispan_y_loc, labels_dict, subsonic_calc_type).get_data(mach_numbers)

    else:
        print("\n****************\nInvalid Entry. Please run script again.\n****************\n")
        quit()

# Print exit statement to verify completion of process
print()
print("Pressure plotting executed successfully. \n")