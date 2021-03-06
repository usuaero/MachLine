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

        # Identify all data including column headers to identify where desired information is located
        Data_complete = np.genfromtxt(file_name, delimiter=",", skip_header=0, dtype=str)
        # Reformate Data_complete to remove multiple quotes
        Data_complete = np.char.replace(Data_complete,'"',"")

        Data = np.genfromtxt(file_name, delimiter=",", skip_header=1, dtype=float)

        # Verify requested point is in the dataset
        if "C_p_inc" not in Data_complete[0,:]:
            print(f'***{"C_p_inc"} is not found in {file_name}. Quitting')
            exit()

        # Iterate over file columns searching for the x and Cp columns
        for i, val in enumerate(Data_complete[0,:]):
            if val == 'Points:0':
                x_loc = i
            if val == "C_p_inc":
                point_loc = i
            
        C_p_inc = Data[:,point_loc]
        x = Data[:,x_loc]

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

    def form_comparison(self, formulations):
        
        # Initialize lists to differentiate plotting markers
        opacity = [.1, .3, .7, 1]


        # Iterate over the percent semispan locations
        for percent_semispan in self.locations:

            # iterate over applicable angles of attack
            for AoA in self.locations[percent_semispan]:
                
                # Create a list conaining the location of data going to be compared
                data_locations = []
                # Iterate over node counts # This is where plotting will happen, to plot the difference of all the nodal results
                for i, Node in enumerate(self.Nodes):
                    data_locations.append([])

                    # Iterate over formulations
                    for l,form in enumerate(formulations):
                        if form =="source-free":
                            form_adj = "source_free"
                            markersize=10
                        else:
                            form_adj = form
                            markersize=5

                        # Determine the data file and create a label for it    
                        data_file =  percent_semispan + "_percent_semispan/" + AoA + "_degrees_AoA/" + percent_semispan + "_percent_semispan_" + Node + "_nodes_results_" + AoA + "_deg_AoA_" + form_adj + "_formulation.csv"
                        data_complete = np.genfromtxt(data_file, delimiter=",", skip_header=0, dtype=str)
                        data_complete = np.char.replace(data_complete,'"',"") # removes additional "" around column headers
                        data = np.genfromtxt(data_file, delimiter=",", skip_header=1, dtype=float)
                        Curve_label = Node + " Nodes, " + form

                        # Verify x and y labels are contained in datafile
                        if ("C_p_inc" not in data_complete[0,:]) or ("arc_length" not in data_complete[0,:]):
                            print(f"***Identified axis labels for formulation comparison and not found in {data_file}")
                            print("Quitting...")
                            exit()

                        # Locate required columns of data in the data_file
                        for k,column in enumerate(data_complete[0,:]):
                            if "C_p_inc" == column:
                                y_loc = k
                            if column == "arc_length":
                                x_loc = k

                        # Pull and reshape data from data file for ease of plotting and organization
                        x_data = data[:,x_loc]
                        y_data = data[:, y_loc]
                        x = np.reshape(x_data, (len(x_data),1))
                        y = np.reshape(y_data, (len(y_data),1))
                        plot_data = np.hstack((x, y))
                        
                        # Assign organized data to storage list
                        data_locations[i].append(plot_data)

                        plt.plot(data_locations[i][l][:,0], data_locations[i][l][:,1], label=Curve_label, color="k", marker="o", markersize=markersize, alpha=opacity[i], fillstyle="none", linestyle="none")
                
                # Format the plot and show
                plt.xlabel("Arc Length")
                plt.ylabel("$C_p$")
                plt.gca().invert_yaxis()
                plt.legend()
                plt.title(AoA + " deg, " + percent_semispan + " percent semispan")
                plt.show()
                # breakpoint()

                        






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

# Proccesses either one or all semispan locations based on input file
if json_vals["plots"]["process_all"]:
    for form in formulation:

        if json_vals["plots"]["display pressure plots"]:
            Swept_Plotting(locations, chord, Nodes, semispan_xy_loc).get_data(form)

    # Compare formulation results
    if json_vals["plots"]["compare formulations"]:
        Swept_Plotting(locations, chord, Nodes, semispan_xy_loc).form_comparison(formulation)

else:
    x = input("Enter percent semispan for analysis results. Options are 4.1, 8.2, 16.3, 24.5, 36.7, 51.0, 65.3, 89.8, 94.9:   ",)
    for form in formulation:
        if x == "51":
            x = "51.0"

        if x in locations.keys():

            Specific_Semispan = {x: locations[x]}
            if json_vals["plots"]["display pressure plots"]:
                Swept_Plotting(Specific_Semispan, chord, Nodes, semispan_xy_loc).get_data(form)

        else:
            print("\n****************\nInvalid Entry. Please run script again.\n****************\n")
            quit()
    
    # Compare formulation results
    if json_vals["plots"]["compare formulations"]:
        Swept_Plotting(Specific_Semispan, chord, Nodes, semispan_xy_loc).form_comparison(formulation)

    # Print exit statement to verify completion of process
    print()
    print("Pressure plotting executed successfully. \n")

