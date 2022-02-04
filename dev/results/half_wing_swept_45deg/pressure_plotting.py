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

            colors = 'g'

        elif Node_Count == '17':

            colors = 'b'

        elif Node_Count == '25':

            colors = 'c'

        elif Node_Count == '45':

            colors = 'r'

        else:
            colors = 'y'
        
        plt.plot((x-LE_xy_loc)/Chord_Length,C_p_inc, label=Plot_label, color = colors)

        return



    def get_data(self):

        # Pull In Results Depending on Location Input
        # In addition to defining leading edge location and respective labels for plot
        for percent_semispan in self.locations:

            # Pull in experimental results
            file_experimental = percent_semispan + "_percent_semispan/" + percent_semispan + "_percent_semispan_Experimental_Results_Wing_A.csv"
        
            Experimental = np.genfromtxt(file_experimental, delimiter=",", skip_header=2, dtype=float)
            
            # Constrain which node counts are used near tip conditions
            if percent_semispan == '89.8':
                self.Nodes = ['17', '25', '45']

            elif percent_semispan == '94.9':
                self.Nodes = ['45']

            # Establish clean list for saving and combining produced plots
            list_of_files =[]
            # Iterate over angles of attack for each semispan
            for j,AoA in enumerate(self.locations[percent_semispan]):


                # Iterate over node counts
                for Node in self.Nodes:   
                    file =  percent_semispan + "_percent_semispan/" + AoA + "_degrees_AoA/" + percent_semispan + "_percent_semispan_" + Node + "_nodes_results_" + AoA + "_deg_AoA.csv"
                    Curve_label = Node + " Nodes"
                    self.Pressure_Plot(file, Node, Curve_label, self.semispan_xy_loc[percent_semispan], self.chord)
                


                #=== Main Plotting Section ===
                #Identify labels and notes to include on plots
                title = r"Pressure Distribution at $\frac{2y}{b} =$ " + str("{:.3f}".format(float(percent_semispan)/100))
               
                #Identify angle of atack note to add to plot title
                y = r'$\alpha$ = '
                Notes = y + AoA + r'$^{\circ}$'

                #Pull in experimental results for comparison at each angle of attack differentiating upper and lower surfaces
                upper_surface_count = Experimental[:,0].size//2 + Experimental[:,0].size%2
                plt.plot(Experimental[:upper_surface_count+1,0], Experimental[:upper_surface_count+1,j+1], ".",color="k", label="Exerimental Upper Surface")
                plt.plot(Experimental[upper_surface_count:,0], Experimental[upper_surface_count:,j+1], "*",color="k", label="Exerimental Lower Surface")

                #Plot the figure containing all curves
                complete_title = title + ", " + Notes
                plt.title(complete_title)
                plt.xlabel('x/c')
                plt.ylabel(r"$C_p$")
                plt.gca().invert_yaxis()
                if percent_semispan == '94.9':
                    plt.legend()
                else:
                    plt.legend(ncol = 2)
                #Save the figure in appropriate location
                
                if json_vals["plots"]["save plots"]:
                    filename = percent_semispan + "_percent_semispan/plots_" + percent_semispan + "_percent_semispan/" + json_vals["plots"]["save plot type"] + "_plots/" + AoA + "degrees_AoA_plot." + json_vals["plots"]["save plot type"]
                    list_of_files.append(filename)
                    plt.savefig(filename)

                plt.show()

            # Combine plots into one pdf for easy review
            if json_vals["plots"]["save plot type"] == "pdf" and json_vals["plots"]["combine pdfs"]:

                merger = PdfFileMerger()
                Combined_name = percent_semispan + '_percent_semispan/plots_' + percent_semispan + '_percent_semispan/pdf_plots/Combined_plots_' + percent_semispan + '_percent_semispan.pdf'
                Combined_name_copy = 'Plot_Summary/' + 'Combined_plots_' + percent_semispan + '_percent_semispan.pdf'

                for name in list_of_files:
                    merger.append(name)

                merger.write(Combined_name)
                merger.write(Combined_name_copy)
                merger.close()





# Main script

inputfile = "Swept_half_wing_conditions_input.json"

json_string = open(inputfile).read()
json_vals = json.loads(json_string)

# Identify values to pass from input file
semispan_xy_loc = {
    '4.1': 2, 
    '8.2': 4,
    '16.3': 7.99,
    '24.5': 12,
    '36.7': 17.98,
    '51.0': 24.99,
    '65.3': 32,
    '89.8': 44,
    '94.9': 46.5
}

chord = json_vals["geometry"]["Uniform Chord Length"]

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

Nodes = ['10', '17', '25', '45']

if json_vals["plots"]["process_all"]:

    Swept_Plotting(locations, chord, Nodes, semispan_xy_loc).get_data()

else:
    x = input("Enter percent semispan for analysis results. Options are 4.1, 8.2, 16.3, 24.5, 36.7, 51.0, 65.3, 89.8, 94.9:   ",)
    if x == "51":
        x = "51.0"

    if x in locations.keys():

        Specific_Semispan = {x: locations[x]}

        Swept_Plotting(Specific_Semispan, chord, Nodes, semispan_xy_loc).get_data()

    else:
        print("\n****************\nInvalid Entry. Please run script again.\n****************\n")
        quit()