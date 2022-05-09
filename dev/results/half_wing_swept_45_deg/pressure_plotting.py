#This file is to take .csv or .txt files of lift distributions and plot
#them allowing comparison to experimental results

import numpy as np
import matplotlib.pyplot as plt
import json
from PyPDF2 import PdfFileMerger

class Swept_Plotting:

    def __init__(self, Geometry_Data, chord_length, Node_Counts, xy_locations):
        self.chord = chord_length

        self.locations = Geometry_Data

        self.semispan_xy_loc = xy_locations

        self.Nodes = Node_Counts
       

    def Pressure_Plot(self, file_name, Node_Count, Plot_label,LE_xy_loc, Chord_Length):

        Data = np.genfromtxt(file_name, delimiter=",", skip_header=1, dtype=float)

        C_p_inc = Data[:,0]
        x = Data[:,1]

        # Ensure the plots are colored consistently

        if Node_Count == '10':
            shape = "o"
            opacity = 0.25
            colors = 'g'

        elif Node_Count == '20':
            shape = "x"
            opacity = 0.5
            colors = 'b'

        elif Node_Count == '40':
            shape = "1"
            opacity = 0.7
            colors = 'c'

        elif Node_Count == '80':
            shape = "2"
            opacity = 1
            colors = 'r'

        else:
            shape = "4"   
            opacity = 0.1
            colors = 'y'
        
        # Black and white lines with varying opacity
        plt.plot((x-LE_xy_loc)/Chord_Length,C_p_inc, label=Plot_label, alpha=opacity, color="k")
        
        # Various colors
        # plt.plot((x-LE_xy_loc)/Chord_Length,C_p_inc, "o", label=Plot_label, color=colors)
        
        # Markers alone without lines for all data
        # plt.plot((x-LE_xy_loc)/Chord_Length,C_p_inc, label=Plot_label, marker=shape, color="k", linestyle="none")

        return



    def get_data(self, formulation):

        # Determine formulation formatting for pdf file purposes
        if formulation == "source-free":
            formulation_adjusted = "source_free"
        else:
            formulation_adjusted = formulation

        # Pull In Results Depending on Location Input
        # In addition to defining leading edge location and respective labels for plot
        for percent_semispan in self.locations:

            # Pull in experimental results
            file_experimental = percent_semispan + "_percent_semispan/" + percent_semispan + "_percent_semispan_Experimental_Results_Wing_A.csv"
        
            Experimental = np.genfromtxt(file_experimental, delimiter=",", skip_header=2, dtype=float)

            # Establish clean list for saving and combining produced plots
            list_of_files =[]
            # Iterate over angles of attack for each semispan
            for j,AoA in enumerate(self.locations[percent_semispan]):


                # Iterate over node counts backwards to allow plot sizing based on highest node count
                for k,Node in enumerate(reversed(self.Nodes)): 

                    file =  percent_semispan + "_percent_semispan/" + AoA + "_degrees_AoA/" + percent_semispan + "_percent_semispan_" + Node + "_nodes_results_" + AoA + "_deg_AoA_" + formulation_adjusted + "_formulation.csv"
                    Curve_label = Node + " Nodes"
                    self.Pressure_Plot(file, Node, Curve_label, self.semispan_xy_loc[percent_semispan], self.chord)
                    # Determine sizing for plot formats based on highest node count used
                    if k == 0:
                        
                        xmin, xmax = plt.xlim()
                        ymin, ymax = plt.ylim()
                
                # Scale plot based on highest node count, and invert the Y-axis
                plt.xlim(xmin, xmax)
                plt.ylim(ymin, ymax)
                plt.gca().invert_yaxis()


                #=== Main Plotting Section ===
                #Identify labels and notes to include on plots
                title = r"Pressure Distribution at $\frac{2y}{b}=$" + str("{:.3f}".format(float(percent_semispan)/100))
               
                #Identify angle of atack note to add to plot title
                y = r'$\alpha$='
                AoA_Notes = y + AoA + r'$^{\circ}$'

                #Pull in experimental results for comparison at each angle of attack differentiating upper and lower surfaces
                upper_surface_count = Experimental[:,0].size//2 + Experimental[:,0].size%2
                plt.plot(Experimental[:upper_surface_count+1,0], Experimental[:upper_surface_count+1,j+1], ".",color="k", label="Exerimental Upper Surface", fillstyle="full")
                plt.plot(Experimental[upper_surface_count:,0], Experimental[upper_surface_count:,j+1], "*",color="k", label="Exerimental Lower Surface", fillstyle="full")

                #Plot the figure containing all curves
                complete_title = title + ", " + AoA_Notes
                Sub_Note = "*" + formulation + " form."
                # plt.title(complete_title)
                xlabel = 'x/c'
                plt.xlabel(xlabel)
                plt.ylabel(r"$C_p$")



                # Identify location of formulation type on plot
                # plt.figtext(0.625, 0.015, Sub_Note)

                # Split legend into columns
                plt.legend(ncol = 2, fontsize=6)
                
                #Save the figure in appropriate location and with the correct size
                
                if json_vals["plots"]["save plots"]:
                    filename = percent_semispan + "_percent_semispan/plots_" + percent_semispan + "_percent_semispan/" + json_vals["plots"]["save plot type"] + "_plots/" + AoA + "degrees_AoA_plot_" + formulation_adjusted + "_formulation." + json_vals["plots"]["save plot type"]
                    list_of_files.append(filename)
                    plt.savefig(filename)
                plt.show()

            # Combine plots into one pdf for easy review
            if json_vals["plots"]["save plot type"] == "pdf" and json_vals["plots"]["combine pdfs"]:

                # Combine files in percent semispan folders
                merger = PdfFileMerger()
                Combined_name = percent_semispan + '_percent_semispan/plots_' + percent_semispan + '_percent_semispan/pdf_plots/Combined_plots_' + percent_semispan + '_percent_semispan_' + formulation_adjusted + '_formulation.pdf'

                for name in list_of_files:
                    merger.append(name)

                merger.write(Combined_name)
                merger.close()

                # Combine files in Plot Summary folder
                Combined_name_copy = 'Plot_Summary/' + 'Combined_plots_' + percent_semispan + '_percent_semispan_' + formulation_adjusted + '_formulation.pdf'
                merger = PdfFileMerger()
                for name in list_of_files:
                    merger.append(name)
                merger.write(Combined_name_copy)
                merger.close()

# Main script
inputfile = "Swept_half_wing_conditions_input.json"
json_string = open(inputfile).read()
json_vals = json.loads(json_string)

# Identify values to pass from input file

locations = json_vals["geometry"]["semispan locations and AoA"]
semispan_xy_loc = json_vals["geometry"]["semispan_xy_loc"]

chord = json_vals["geometry"]["uniform chord length"]
Nodes = json_vals["geometry"]["nodes"]
formulation = json_vals["solver"]["formulation"]

# Iterate over solver formulations from input file

    # Proccesses either one or all semispan locations based on input file
if json_vals["plots"]["process_all"]:
    for form in formulation:

        Swept_Plotting(locations, chord, Nodes, semispan_xy_loc).get_data(form)

else:
    x = input("Enter percent semispan for analysis results. Options are 4.1, 8.2, 16.3, 24.5, 36.7, 51.0, 65.3, 89.8, 94.9:   ",)
    for form in formulation:
        if x == "51":
            x = "51.0"

        if x in locations.keys():

            Specific_Semispan = {x: locations[x]}

            Swept_Plotting(Specific_Semispan, chord, Nodes, semispan_xy_loc).get_data(form)

        else:
            print("\n****************\nInvalid Entry. Please run script again.\n****************\n")
            quit()

    # Print exit statement to verify completion of process
    print()
    print("Pressure plotting executed successfully. \n")