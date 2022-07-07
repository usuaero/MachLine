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

        # Iterate over file columns searching for specific columns
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
        
        #===========Adjustable plotting options left for user's convenience==========#

        # Black and white lines with varying opacity
        plt.plot((x-LE_xy_loc)/Chord_Length,C_p_inc, label=Plot_label, alpha=opacity, color="k")
        
        # Various colored lines
        # plt.plot((x-LE_xy_loc)/Chord_Length,C_p_inc, "o", label=Plot_label, color=colors)
        
        # Markers only
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

                # Initialize figure size
                plt.figure(figsize=(4,3))
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
                title = "Pressure Distribution at $\frac{2y}{b}=$" + str("{:.3f}".format(float(percent_semispan)/100))
               
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
                plt.ylabel("$C_p$")



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

    def form_comparison(self, formulations, adjust_xscale, adjust_yscale):
        
        # Initialize lists to differentiate plotting markers
        # opacity = [.1, .3, .7, 1]
        size = [15, 10, 7, 3]
        print()

        # Iterate over the percent semispan locations
        for percent_semispan in self.locations:

            # Create a list to store figure names for combining pdfs
            figure_list = []

            # iterate over applicable angles of attack
            for AoA in self.locations[percent_semispan]:
                
                # Create a list conaining the location of data going to be compared and set size of figure to be produced
                data_locations = []
                plt.figure(figsize=(6.4,4.8))
                # Iterate over node counts # This is where plotting will happen, to plot the difference of all the nodal results
                for i, Node in enumerate(self.Nodes):
                    data_locations.append([])

                    force_comparison = []
                    # Iterate over formulations
                    for l,form in enumerate(formulations):
                        if form =="source-free":
                            form_adj = "source_free"
                            markersize=10
                        else:
                            form_adj = form
                            markersize=5

                        force_comparison.append([])
                        # Determine the data file and create a label for it    
                        data_file =  percent_semispan + "_percent_semispan/" + AoA + "_degrees_AoA/" + percent_semispan + "_percent_semispan_" + Node + "_nodes_results_" + AoA + "_deg_AoA_" + form_adj + "_formulation.csv"
                        data_complete = np.genfromtxt(data_file, delimiter=",", skip_header=0, dtype=str)
                        data_complete = np.char.replace(data_complete,'"',"") # removes additional "" around column headers
                        data = np.genfromtxt(data_file, delimiter=",", skip_header=1, dtype=float)

                        # Verify x and y labels are contained in datafile
                        if ("C_p_inc" not in data_complete[0,:]) or ("arc_length" not in data_complete[0,:]):
                            print(f"***Identified axis labels for formulation comparison and not found in {data_file}")
                            print("Quitting...")
                            exit()

                        # Locate required columns of data in the data_file
                        for k,column in enumerate(data_complete[0,:]):
                            if column == "C_p_inc":
                                y_loc = k
                            if column == "arc_length":
                                x_loc = k
                            if column == "dC_f:0":
                                x_force_loc = k
                            if column == "dC_f:1":
                                y_force_loc = k
                            if column == "dC_f:2":
                                z_force_loc = k

                        # Add check to throw errors
                        if (x_force_loc == y_force_loc) or (x_force_loc == z_force_loc) or (y_force_loc == z_force_loc):
                            print("    Error. The identified locations of force columns in data file are conflicting. Quitting...")
                            exit()

                        # Pull and reshape data from data file for ease of plotting and organization
                        x_data = data[:,x_loc]
                        y_data = data[:, y_loc]
                        x = np.reshape(x_data, (len(x_data),1))
                        y = np.reshape(y_data, (len(y_data),1))
                        plot_data = np.hstack((x, y))

                        # Assign organized data to storage lists for each formulation
                        data_locations[i].append(plot_data)
                        force_comparison[l].append(sum(data[:,x_force_loc])) # Append Cx
                        force_comparison[l].append(sum(data[:,y_force_loc])) # Append Cy
                        force_comparison[l].append(sum(data[:,z_force_loc])) # Append Cz
                        
                        # Compute the difference between formulation solutions
                        if l > 0:
                            # Verify that all x axis values are the same
                            count = 0
                            pos_difference = []
                            neg_difference = []
                            x_axis = np.array(data_locations[i][l][:,0]).tolist()
                            pos_x_axis = []
                            neg_x_axis = []
                
                            # Iterate over values in the x axis
                            for r, val in enumerate(x_axis):

                                # Prevent errors in alignment and comparison of formulation locations
                                if val != data_locations[i][0][r,0]:
                                    
                                    # Ensure error message is only thrown once
                                    count += 1
                                    if count > 1:
                                        print(f"X axis values are not lining up right for the {form} formulation. Results may be invalid.")  
                                
                                # Calculate difference between formulation results
                                diff = data_locations[i][l][r,1] - data_locations[i][0][r,1]
                                
                                if diff >=0.0:
                                    pos_difference.append(diff)
                                    pos_x_axis.append(x_axis[r])
                                else:
                                    neg_difference.append(abs(diff))
                                    neg_x_axis.append(x_axis[r])

                            # Label curve and plot difference    
                            Curve_label = f"{Node} nodes"
                            plt.plot(pos_x_axis, pos_difference, marker="o", markersize=size[i], label=Curve_label, linestyle="none", color="k")
                            
                            # Plot negative difference values a different color
                            if len(neg_difference) != 0:
                                plt.plot(neg_x_axis, neg_difference, marker="o", markersize=size[i], linestyle="none", color="k", alpha=0.3)

                    # Print force differences if verbose is selected
                    if json_vals["verbose"]:
                        # Calculate force difference between formulations at each node count
                        force_diff = np.array(force_comparison[1]) - np.array(force_comparison[0])

                        # Print force results to terminal to reduce clutter on plots
                        if i == 0: print(f"Force differences for {percent_semispan} percent semispan and {AoA} deg AoA")
                        print(f'{Node} nodes:')
                        print(f"    Cx({formulation[1]} - {formulation[0]}):   ", force_diff[0])
                        print(f"    Cy({formulation[1]} - {formulation[0]}):   ", force_diff[1])
                        print(f"    Cz({formulation[1]} - {formulation[0]}):   ", force_diff[2])

                # Format the plot
                if adjust_xscale != "auto":
                    plt.xlim(-adjust_xscale, adjust_xscale)
                if adjust_yscale != "auto":
                    plt.ylim(-adjust_yscale, adjust_yscale)
                plt.xlabel("Arc Length")
                plt.ylabel(f"$|\Delta C_p|$ ({formulations[1]} - {formulations[0]})")
                # plt.gca().invert_yaxis()
                plt.yscale("log")
                plt.legend()

                # Include title if verbose is selected
                if json_vals["verbose"]:
                    title = "Formulation Comparison at " + percent_semispan + " percent semispan"
                    angle = r'$\alpha$=' + AoA + r'$^{\circ}$'

                    plt.figtext(0.825, 0.025, angle)
                    plt.title(title)
                    print()

                # Save the comparison plots in the appropriate locations
                if json_vals["plots"]["save plots"]:
                    comparison_plot_name = percent_semispan + "_percent_semispan/plots_" + percent_semispan + "_percent_semispan/" + json_vals["plots"]["save plot type"] + "_plots/Formulation comparison plot " + AoA + " degrees_AoA." + json_vals["plots"]["save plot type"]
                    figure_list.append(comparison_plot_name)
                    plt.savefig(comparison_plot_name)
                plt.show()
            
            if json_vals["plots"]["combine pdfs"] and json_vals["plots"]["save plot type"] == "pdf":
                # Combine files in plot summary folder
                merger = PdfFileMerger()
                Combined_name = "Plot_Summary/Formulation_comparison_plots_" + percent_semispan + "_percent_semispan.pdf"

                for name in figure_list:
                    merger.append(name)

                merger.write(Combined_name)
                merger.close()

# =======================================================================================================================================================================================================================================================================================
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
        # Verify more than one formulation has been selected
        if len(formulation) < 2:
            print("Enter at least two formulations for comparison to be possible. Quitting...")
            exit()
        xscale = json_vals["plots"]["x axis scale"]
        yscale = json_vals["plots"]["y axis scale"]
        Swept_Plotting(locations, chord, Nodes, semispan_xy_loc).form_comparison(formulation, xscale, yscale)

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
        # Verify more than one formulation has been selected
        if len(formulation) < 2:
            print("Enter at least two formulations for comparison to be possible. Quitting...")
            exit()

        xscale = json_vals["plots"]["x axis scale"]
        yscale = json_vals["plots"]["y axis scale"]
        Swept_Plotting(Specific_Semispan, chord, Nodes, semispan_xy_loc).form_comparison(formulation, xscale, yscale)

# Print exit statement to verify completion of process 
print("Pressure plotting executed successfully. \n")

