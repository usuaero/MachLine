#This file is to take .csv or .txt files of lift distributions and plot
#them allowing comparison to experimental results

import numpy as np
import matplotlib.pyplot as plt
import json
#from PyPDF2 import PdfFileMerger

class Swept_Plotting:

    def __init__(self, LE_Location):
        self.LE_Location = LE_Location

    def Pressure_Plot(self, file, Plot_label, LE_loc, Chord_Length):
        Data = np.genfromtxt(file, delimiter=",", skip_header=1, dtype=float)

        C_p_inc = Data[:,0]
        x = Data[:,1]

        # Ensure the plots are colored consistently
        if LE_loc == '89.8':
                if file == self.file17_first_AoA or file == self.file17_second_AoA or file == self.file17_third_AoA:

                    colors = 'b'

                elif file == self.file25_first_AoA or file == self.file25_second_AoA or file == self.file25_third_AoA:

                    colors = 'c'

                elif file == self.file45_first_AoA or file == self.file45_second_AoA or file == self.file45_third_AoA:

                    colors = 'r'

                else:
                    colors = 'y'
                
        elif LE_loc == '94.9':
                if file == self.file45_first_AoA or file == self.file45_second_AoA or file == self.file45_third_AoA or file == self.file45_fourth_AoA or file == self.file45_fifth_AoA or file == self.file45_sixth_AoA:

                    colors = 'r'

                else:
                    colors = 'y'

        else:
            if file == self.file10_first_AoA or file == self.file10_second_AoA or file == self.file10_third_AoA or file == self.file10_fourth_AoA or file == self.file10_fifth_AoA or file == self.file10_sixth_AoA:

                colors = 'g'

            elif file == self.file17_first_AoA or file == self.file17_second_AoA or file == self.file17_third_AoA or file == self.file17_fourth_AoA or file == self.file17_fifth_AoA or file == self.file17_sixth_AoA:

                colors = 'b'

            elif file == self.file25_first_AoA or file == self.file25_second_AoA or file == self.file25_third_AoA or file == self.file25_fourth_AoA or file == self.file25_fifth_AoA or file == self.file25_sixth_AoA:

                colors = 'c'

            elif file == self.file45_first_AoA or file == self.file45_second_AoA or file == self.file45_third_AoA or file == self.file45_fourth_AoA or file == self.file45_fifth_AoA or file == self.file45_sixth_AoA:

                colors = 'r'

            else:
                colors = 'y'
        
        plt.plot((x-self.x_LE)/Chord_Length,C_p_inc, label=Plot_label, color = colors)

        return


   # LE_Location = input("Enter percent semispan for analysis results. Options are 4.1, 8.2, 16.3, 24.5, 36.7, 51.0, 65.3, 89.8, 94.9:   ",)

    def get_data(self):
        chord = 20
        #Pull In Results Depending on Location Input
        #In addition to defining leading edge location and respective labels for plot
        
        if self.LE_Location == '4.1':
            #File locations at an angle of attack of 0 degrees   
            self.file10_first_AoA = "4.1_percent_semispan/0_degrees_AoA/4.1_percent_semispan_10_nodes_results_0_deg_AoA.csv"
            self.file17_first_AoA = "4.1_percent_semispan/0_degrees_AoA/4.1_percent_semispan_17_nodes_results_0_deg_AoA.csv"
            self.file25_first_AoA = "4.1_percent_semispan/0_degrees_AoA/4.1_percent_semispan_25_nodes_results_0_deg_AoA.csv"
            self.file45_first_AoA = "4.1_percent_semispan/0_degrees_AoA/4.1_percent_semispan_45_nodes_results_0_deg_AoA.csv"
            
            #File locations at an angle of attack of 2.1 degrees
            self.file10_second_AoA = "4.1_percent_semispan/2.1_degrees_AoA/4.1_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
            self.file17_second_AoA = "4.1_percent_semispan/2.1_degrees_AoA/4.1_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
            self.file25_second_AoA = "4.1_percent_semispan/2.1_degrees_AoA/4.1_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
            self.file45_second_AoA = "4.1_percent_semispan/2.1_degrees_AoA/4.1_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

            #File locations at an angle of attack of 4.2 degrees
            self.file10_third_AoA = "4.1_percent_semispan/4.2_degrees_AoA/4.1_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
            self.file17_third_AoA = "4.1_percent_semispan/4.2_degrees_AoA/4.1_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
            self.file25_third_AoA = "4.1_percent_semispan/4.2_degrees_AoA/4.1_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
            self.file45_third_AoA = "4.1_percent_semispan/4.2_degrees_AoA/4.1_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

            #File locations at an angle of attack of 6.2 degrees
            self.file10_fourth_AoA = "4.1_percent_semispan/6.2_degrees_AoA/4.1_percent_semispan_10_nodes_results_6.2_deg_AoA.csv"
            self.file17_fourth_AoA = "4.1_percent_semispan/6.2_degrees_AoA/4.1_percent_semispan_17_nodes_results_6.2_deg_AoA.csv"
            self.file25_fourth_AoA = "4.1_percent_semispan/6.2_degrees_AoA/4.1_percent_semispan_25_nodes_results_6.2_deg_AoA.csv"
            self.file45_fourth_AoA = "4.1_percent_semispan/6.2_degrees_AoA/4.1_percent_semispan_45_nodes_results_6.2_deg_AoA.csv"

            #File locations at an angle of attack of 8.3 degrees
            self.file10_fifth_AoA = "4.1_percent_semispan/8.3_degrees_AoA/4.1_percent_semispan_10_nodes_results_8.3_deg_AoA.csv"
            self.file17_fifth_AoA = "4.1_percent_semispan/8.3_degrees_AoA/4.1_percent_semispan_17_nodes_results_8.3_deg_AoA.csv"
            self.file25_fifth_AoA = "4.1_percent_semispan/8.3_degrees_AoA/4.1_percent_semispan_25_nodes_results_8.3_deg_AoA.csv"
            self.file45_fifth_AoA = "4.1_percent_semispan/8.3_degrees_AoA/4.1_percent_semispan_45_nodes_results_8.3_deg_AoA.csv"

            #File locations at an angle of attack of 10.4 degrees
            self.file10_sixth_AoA = "4.1_percent_semispan/10.4_degrees_AoA/4.1_percent_semispan_10_nodes_results_10.4_deg_AoA.csv"
            self.file17_sixth_AoA = "4.1_percent_semispan/10.4_degrees_AoA/4.1_percent_semispan_17_nodes_results_10.4_deg_AoA.csv"
            self.file25_sixth_AoA = "4.1_percent_semispan/10.4_degrees_AoA/4.1_percent_semispan_25_nodes_results_10.4_deg_AoA.csv"
            self.file45_sixth_AoA = "4.1_percent_semispan/10.4_degrees_AoA/4.1_percent_semispan_45_nodes_results_10.4_deg_AoA.csv"

            #Location of the experimental results
            file_experimental = "4.1_percent_semispan/4.1_percent_semispan_Experimental_Results_Wing_A.csv"
            
            self.x_LE = 2
            AoA = ['0', '2.1', '4.2', '6.2', '8.3', '10.4']
            
            #Identify labels and notes to include on plots
            labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
            title = r"Pressure Coefficient at $\frac{2y}{b} = 0.041$"
                

        elif self.LE_Location == '8.2':
            #File locations at an angle of attack of 0 degrees
            self.file10_first_AoA = "8.2_percent_semispan/0_degrees_AoA/8.2_percent_semispan_10_nodes_results_0_deg_AoA.csv"
            self.file17_first_AoA = "8.2_percent_semispan/0_degrees_AoA/8.2_percent_semispan_17_nodes_results_0_deg_AoA.csv"
            self.file25_first_AoA = "8.2_percent_semispan/0_degrees_AoA/8.2_percent_semispan_25_nodes_results_0_deg_AoA.csv"
            self.file45_first_AoA = "8.2_percent_semispan/0_degrees_AoA/8.2_percent_semispan_45_nodes_results_0_deg_AoA.csv"

            #File locations at an angle of attack of 2.1 degrees
            self.file10_second_AoA = "8.2_percent_semispan/2.1_degrees_AoA/8.2_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
            self.file17_second_AoA = "8.2_percent_semispan/2.1_degrees_AoA/8.2_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
            self.file25_second_AoA = "8.2_percent_semispan/2.1_degrees_AoA/8.2_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
            self.file45_second_AoA = "8.2_percent_semispan/2.1_degrees_AoA/8.2_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

            #File locations at an angle of attack of 4.2 degrees
            self.file10_third_AoA = "8.2_percent_semispan/4.2_degrees_AoA/8.2_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
            self.file17_third_AoA = "8.2_percent_semispan/4.2_degrees_AoA/8.2_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
            self.file25_third_AoA = "8.2_percent_semispan/4.2_degrees_AoA/8.2_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
            self.file45_third_AoA = "8.2_percent_semispan/4.2_degrees_AoA/8.2_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

            #File locations at an angle of attack of 6.2 degrees
            self.file10_fourth_AoA = "8.2_percent_semispan/6.2_degrees_AoA/8.2_percent_semispan_10_nodes_results_6.2_deg_AoA.csv"
            self.file17_fourth_AoA = "8.2_percent_semispan/6.2_degrees_AoA/8.2_percent_semispan_17_nodes_results_6.2_deg_AoA.csv"
            self.file25_fourth_AoA = "8.2_percent_semispan/6.2_degrees_AoA/8.2_percent_semispan_25_nodes_results_6.2_deg_AoA.csv"
            self.file45_fourth_AoA = "8.2_percent_semispan/6.2_degrees_AoA/8.2_percent_semispan_45_nodes_results_6.2_deg_AoA.csv"

            #File locations at an angle of attack of 8.3 degrees
            self.file10_fifth_AoA = "8.2_percent_semispan/8.3_degrees_AoA/8.2_percent_semispan_10_nodes_results_8.3_deg_AoA.csv"
            self.file17_fifth_AoA = "8.2_percent_semispan/8.3_degrees_AoA/8.2_percent_semispan_17_nodes_results_8.3_deg_AoA.csv"
            self.file25_fifth_AoA = "8.2_percent_semispan/8.3_degrees_AoA/8.2_percent_semispan_25_nodes_results_8.3_deg_AoA.csv"
            self.file45_fifth_AoA = "8.2_percent_semispan/8.3_degrees_AoA/8.2_percent_semispan_45_nodes_results_8.3_deg_AoA.csv"

            #File locations at an angle of attack of 10.4 degrees
            self.file10_sixth_AoA = "8.2_percent_semispan/10.4_degrees_AoA/8.2_percent_semispan_10_nodes_results_10.4_deg_AoA.csv"
            self.file17_sixth_AoA = "8.2_percent_semispan/10.4_degrees_AoA/8.2_percent_semispan_17_nodes_results_10.4_deg_AoA.csv"
            self.file25_sixth_AoA = "8.2_percent_semispan/10.4_degrees_AoA/8.2_percent_semispan_25_nodes_results_10.4_deg_AoA.csv"
            self.file45_sixth_AoA = "8.2_percent_semispan/10.4_degrees_AoA/8.2_percent_semispan_45_nodes_results_10.4_deg_AoA.csv"

            #Location of the experimental results
            file_experimental = "8.2_percent_semispan/8.2_percent_semispan_Experimental_Results_Wing_A.csv"

            self.x_LE = 4
            AoA = ['0', '2.1', '4.2', '6.2', '8.3', '10.4']
            
            #Identify labels and title to include on plots
            labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
            title = r"Pressure Coefficient at $\frac{2y}{b} = 0.082$"


        elif self.LE_Location == '16.3':
            #File locations at an angle of attack of 0 degrees
            self.file10_first_AoA = "16.3_percent_semispan/0_degrees_AoA/16.3_percent_semispan_10_nodes_results_0_deg_AoA.csv"
            self.file17_first_AoA = "16.3_percent_semispan/0_degrees_AoA/16.3_percent_semispan_17_nodes_results_0_deg_AoA.csv"
            self.file25_first_AoA = "16.3_percent_semispan/0_degrees_AoA/16.3_percent_semispan_25_nodes_results_0_deg_AoA.csv"
            self.file45_first_AoA = "16.3_percent_semispan/0_degrees_AoA/16.3_percent_semispan_45_nodes_results_0_deg_AoA.csv"
            
            #File locations at an angle of attack of 2.1 degrees
            self.file10_second_AoA = "16.3_percent_semispan/2.1_degrees_AoA/16.3_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
            self.file17_second_AoA = "16.3_percent_semispan/2.1_degrees_AoA/16.3_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
            self.file25_second_AoA = "16.3_percent_semispan/2.1_degrees_AoA/16.3_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
            self.file45_second_AoA = "16.3_percent_semispan/2.1_degrees_AoA/16.3_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

            #File locations at an angle of attack of 4.2 degrees
            self.file10_third_AoA = "16.3_percent_semispan/4.2_degrees_AoA/16.3_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
            self.file17_third_AoA = "16.3_percent_semispan/4.2_degrees_AoA/16.3_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
            self.file25_third_AoA = "16.3_percent_semispan/4.2_degrees_AoA/16.3_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
            self.file45_third_AoA = "16.3_percent_semispan/4.2_degrees_AoA/16.3_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

            #File locations at an angle of attack of 6.2 degrees
            self.file10_fourth_AoA = "16.3_percent_semispan/6.2_degrees_AoA/16.3_percent_semispan_10_nodes_results_6.2_deg_AoA.csv"
            self.file17_fourth_AoA = "16.3_percent_semispan/6.2_degrees_AoA/16.3_percent_semispan_17_nodes_results_6.2_deg_AoA.csv"
            self.file25_fourth_AoA = "16.3_percent_semispan/6.2_degrees_AoA/16.3_percent_semispan_25_nodes_results_6.2_deg_AoA.csv"
            self.file45_fourth_AoA = "16.3_percent_semispan/6.2_degrees_AoA/16.3_percent_semispan_45_nodes_results_6.2_deg_AoA.csv"

            #File locations at an angle of attack of 8.3 degrees
            self.file10_fifth_AoA = "16.3_percent_semispan/8.3_degrees_AoA/16.3_percent_semispan_10_nodes_results_8.3_deg_AoA.csv"
            self.file17_fifth_AoA = "16.3_percent_semispan/8.3_degrees_AoA/16.3_percent_semispan_17_nodes_results_8.3_deg_AoA.csv"
            self.file25_fifth_AoA = "16.3_percent_semispan/8.3_degrees_AoA/16.3_percent_semispan_25_nodes_results_8.3_deg_AoA.csv"
            self.file45_fifth_AoA = "16.3_percent_semispan/8.3_degrees_AoA/16.3_percent_semispan_45_nodes_results_8.3_deg_AoA.csv"

            #File locations at an angle of attack of 10.4 degrees
            self.file10_sixth_AoA = "16.3_percent_semispan/10.4_degrees_AoA/16.3_percent_semispan_10_nodes_results_10.4_deg_AoA.csv"
            self.file17_sixth_AoA = "16.3_percent_semispan/10.4_degrees_AoA/16.3_percent_semispan_17_nodes_results_10.4_deg_AoA.csv"
            self.file25_sixth_AoA = "16.3_percent_semispan/10.4_degrees_AoA/16.3_percent_semispan_25_nodes_results_10.4_deg_AoA.csv"
            self.file45_sixth_AoA = "16.3_percent_semispan/10.4_degrees_AoA/16.3_percent_semispan_45_nodes_results_10.4_deg_AoA.csv"

            #Location of the experimental results
            file_experimental = "16.3_percent_semispan/16.3_percent_semispan_Experimental_Results_Wing_A.csv"
            
            self.x_LE = 7.99
            AoA = ['0', '2.1', '4.2', '6.2', '8.3', '10.4']
            
            #Identify labels and title to include on plots
            labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
            title = r"Pressure Coefficient at $\frac{2y}{b} = 0.163$"


        elif self.LE_Location == '24.5':
            #File locations at an angle of attack of 0 degrees
            self.file10_first_AoA = "24.5_percent_semispan/0_degrees_AoA/24.5_percent_semispan_10_nodes_results_0_deg_AoA.csv"
            self.file17_first_AoA = "24.5_percent_semispan/0_degrees_AoA/24.5_percent_semispan_17_nodes_results_0_deg_AoA.csv"
            self.file25_first_AoA = "24.5_percent_semispan/0_degrees_AoA/24.5_percent_semispan_25_nodes_results_0_deg_AoA.csv"
            self.file45_first_AoA = "24.5_percent_semispan/0_degrees_AoA/24.5_percent_semispan_45_nodes_results_0_deg_AoA.csv"

            #File locations at an angle of attack of 2.1 degrees
            self.file10_second_AoA = "24.5_percent_semispan/2.1_degrees_AoA/24.5_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
            self.file17_second_AoA = "24.5_percent_semispan/2.1_degrees_AoA/24.5_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
            self.file25_second_AoA = "24.5_percent_semispan/2.1_degrees_AoA/24.5_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
            self.file45_second_AoA = "24.5_percent_semispan/2.1_degrees_AoA/24.5_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

            #File locations at an angle of attack of 4.2 degrees
            self.file10_third_AoA = "24.5_percent_semispan/4.2_degrees_AoA/24.5_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
            self.file17_third_AoA = "24.5_percent_semispan/4.2_degrees_AoA/24.5_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
            self.file25_third_AoA = "24.5_percent_semispan/4.2_degrees_AoA/24.5_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
            self.file45_third_AoA = "24.5_percent_semispan/4.2_degrees_AoA/24.5_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

            #File locations at an angle of attack of 6.2 degrees
            self.file10_fourth_AoA = "24.5_percent_semispan/6.2_degrees_AoA/24.5_percent_semispan_10_nodes_results_6.2_deg_AoA.csv"
            self.file17_fourth_AoA = "24.5_percent_semispan/6.2_degrees_AoA/24.5_percent_semispan_17_nodes_results_6.2_deg_AoA.csv"
            self.file25_fourth_AoA = "24.5_percent_semispan/6.2_degrees_AoA/24.5_percent_semispan_25_nodes_results_6.2_deg_AoA.csv"
            self.file45_fourth_AoA = "24.5_percent_semispan/6.2_degrees_AoA/24.5_percent_semispan_45_nodes_results_6.2_deg_AoA.csv"

            #File locations at an angle of attack of 8.3 degrees
            self.file10_fifth_AoA = "24.5_percent_semispan/8.3_degrees_AoA/24.5_percent_semispan_10_nodes_results_8.3_deg_AoA.csv"
            self.file17_fifth_AoA = "24.5_percent_semispan/8.3_degrees_AoA/24.5_percent_semispan_17_nodes_results_8.3_deg_AoA.csv"
            self.file25_fifth_AoA = "24.5_percent_semispan/8.3_degrees_AoA/24.5_percent_semispan_25_nodes_results_8.3_deg_AoA.csv"
            self.file45_fifth_AoA = "24.5_percent_semispan/8.3_degrees_AoA/24.5_percent_semispan_45_nodes_results_8.3_deg_AoA.csv"

            #File locations at an angle of attack of 10.4 degrees
            self.file10_sixth_AoA = "24.5_percent_semispan/10.4_degrees_AoA/24.5_percent_semispan_10_nodes_results_10.4_deg_AoA.csv"
            self.file17_sixth_AoA = "24.5_percent_semispan/10.4_degrees_AoA/24.5_percent_semispan_17_nodes_results_10.4_deg_AoA.csv"
            self.file25_sixth_AoA = "24.5_percent_semispan/10.4_degrees_AoA/24.5_percent_semispan_25_nodes_results_10.4_deg_AoA.csv"
            self.file45_sixth_AoA = "24.5_percent_semispan/10.4_degrees_AoA/24.5_percent_semispan_45_nodes_results_10.4_deg_AoA.csv"

            #Location of the experimental results
            file_experimental = "24.5_percent_semispan/24.5_percent_semispan_Experimental_Results_Wing_A.csv"

            self.x_LE = 12
            AoA = ['0', '2.1', '4.2', '6.2', '8.3', '10.4']

            labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
            title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.245$")


        elif self.LE_Location == '36.7':
            #File locations at an angle of attack of 0 degrees
            self.file10_first_AoA = "36.7_percent_semispan/0_degrees_AoA/36.7_percent_semispan_10_nodes_results_0_deg_AoA.csv"
            self.file17_first_AoA = "36.7_percent_semispan/0_degrees_AoA/36.7_percent_semispan_17_nodes_results_0_deg_AoA.csv"
            self.file25_first_AoA = "36.7_percent_semispan/0_degrees_AoA/36.7_percent_semispan_25_nodes_results_0_deg_AoA.csv"
            self.file45_first_AoA = "36.7_percent_semispan/0_degrees_AoA/36.7_percent_semispan_45_nodes_results_0_deg_AoA.csv"

            #File locations at an angle of attack of 2.1 degrees
            self.file10_second_AoA = "36.7_percent_semispan/2.1_degrees_AoA/36.7_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
            self.file17_second_AoA = "36.7_percent_semispan/2.1_degrees_AoA/36.7_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
            self.file25_second_AoA = "36.7_percent_semispan/2.1_degrees_AoA/36.7_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
            self.file45_second_AoA = "36.7_percent_semispan/2.1_degrees_AoA/36.7_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

            #File locations at an angle of attack of 4.2 degrees
            self.file10_third_AoA = "36.7_percent_semispan/4.2_degrees_AoA/36.7_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
            self.file17_third_AoA = "36.7_percent_semispan/4.2_degrees_AoA/36.7_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
            self.file25_third_AoA = "36.7_percent_semispan/4.2_degrees_AoA/36.7_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
            self.file45_third_AoA = "36.7_percent_semispan/4.2_degrees_AoA/36.7_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

            #File locations at an angle of attack of 6.3 degrees
            self.file10_fourth_AoA = "36.7_percent_semispan/6.3_degrees_AoA/36.7_percent_semispan_10_nodes_results_6.3_deg_AoA.csv"
            self.file17_fourth_AoA = "36.7_percent_semispan/6.3_degrees_AoA/36.7_percent_semispan_17_nodes_results_6.3_deg_AoA.csv"
            self.file25_fourth_AoA = "36.7_percent_semispan/6.3_degrees_AoA/36.7_percent_semispan_25_nodes_results_6.3_deg_AoA.csv"
            self.file45_fourth_AoA = "36.7_percent_semispan/6.3_degrees_AoA/36.7_percent_semispan_45_nodes_results_6.3_deg_AoA.csv"

            #File locations at an angle of attack of 8.4 degrees
            self.file10_fifth_AoA = "36.7_percent_semispan/8.4_degrees_AoA/36.7_percent_semispan_10_nodes_results_8.4_deg_AoA.csv"
            self.file17_fifth_AoA = "36.7_percent_semispan/8.4_degrees_AoA/36.7_percent_semispan_17_nodes_results_8.4_deg_AoA.csv"
            self.file25_fifth_AoA = "36.7_percent_semispan/8.4_degrees_AoA/36.7_percent_semispan_25_nodes_results_8.4_deg_AoA.csv"
            self.file45_fifth_AoA = "36.7_percent_semispan/8.4_degrees_AoA/36.7_percent_semispan_45_nodes_results_8.4_deg_AoA.csv"

            #File locations at an angle of attack of 10.5 degrees
            self.file10_sixth_AoA = "36.7_percent_semispan/10.5_degrees_AoA/36.7_percent_semispan_10_nodes_results_10.5_deg_AoA.csv"
            self.file17_sixth_AoA = "36.7_percent_semispan/10.5_degrees_AoA/36.7_percent_semispan_17_nodes_results_10.5_deg_AoA.csv"
            self.file25_sixth_AoA = "36.7_percent_semispan/10.5_degrees_AoA/36.7_percent_semispan_25_nodes_results_10.5_deg_AoA.csv"
            self.file45_sixth_AoA = "36.7_percent_semispan/10.5_degrees_AoA/36.7_percent_semispan_45_nodes_results_10.5_deg_AoA.csv"

            #Location of the experimental results
            file_experimental = "36.7_percent_semispan/36.7_percent_semispan_Experimental_Results_Wing_A.csv"

            self.x_LE = 17.98
            AoA = ['0', '2.1', '4.2', '6.3', '8.4', '10.5']

            labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
            title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.367$")


        elif self.LE_Location == '51.0' or self.LE_Location == '51':
            self.LE_Location = "51.0"

            #File locations at an angle of attack of 0 degrees
            self.file10_first_AoA = "51.0_percent_semispan/0_degrees_AoA/51.0_percent_semispan_10_nodes_results_0_deg_AoA.csv"
            self.file17_first_AoA = "51.0_percent_semispan/0_degrees_AoA/51.0_percent_semispan_17_nodes_results_0_deg_AoA.csv"
            self.file25_first_AoA = "51.0_percent_semispan/0_degrees_AoA/51.0_percent_semispan_25_nodes_results_0_deg_AoA.csv"
            self.file45_first_AoA = "51.0_percent_semispan/0_degrees_AoA/51.0_percent_semispan_45_nodes_results_0_deg_AoA.csv"

            #File locations at an angle of attack of 2.1 degrees
            self.file10_second_AoA = "51.0_percent_semispan/2.1_degrees_AoA/51.0_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
            self.file17_second_AoA = "51.0_percent_semispan/2.1_degrees_AoA/51.0_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
            self.file25_second_AoA = "51.0_percent_semispan/2.1_degrees_AoA/51.0_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
            self.file45_second_AoA = "51.0_percent_semispan/2.1_degrees_AoA/51.0_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

            #File locations at an angle of attack of 4.2 degrees
            self.file10_third_AoA = "51.0_percent_semispan/4.2_degrees_AoA/51.0_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
            self.file17_third_AoA = "51.0_percent_semispan/4.2_degrees_AoA/51.0_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
            self.file25_third_AoA = "51.0_percent_semispan/4.2_degrees_AoA/51.0_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
            self.file45_third_AoA = "51.0_percent_semispan/4.2_degrees_AoA/51.0_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

            #File locations at an angle of attack of 6.3 degrees
            self.file10_fourth_AoA = "51.0_percent_semispan/6.3_degrees_AoA/51.0_percent_semispan_10_nodes_results_6.3_deg_AoA.csv"
            self.file17_fourth_AoA = "51.0_percent_semispan/6.3_degrees_AoA/51.0_percent_semispan_17_nodes_results_6.3_deg_AoA.csv"
            self.file25_fourth_AoA = "51.0_percent_semispan/6.3_degrees_AoA/51.0_percent_semispan_25_nodes_results_6.3_deg_AoA.csv"
            self.file45_fourth_AoA = "51.0_percent_semispan/6.3_degrees_AoA/51.0_percent_semispan_45_nodes_results_6.3_deg_AoA.csv"

            #File locations at an angle of attack of 8.4 degrees
            self.file10_fifth_AoA = "51.0_percent_semispan/8.4_degrees_AoA/51.0_percent_semispan_10_nodes_results_8.4_deg_AoA.csv"
            self.file17_fifth_AoA = "51.0_percent_semispan/8.4_degrees_AoA/51.0_percent_semispan_17_nodes_results_8.4_deg_AoA.csv"
            self.file25_fifth_AoA = "51.0_percent_semispan/8.4_degrees_AoA/51.0_percent_semispan_25_nodes_results_8.4_deg_AoA.csv"
            self.file45_fifth_AoA = "51.0_percent_semispan/8.4_degrees_AoA/51.0_percent_semispan_45_nodes_results_8.4_deg_AoA.csv"

            #File locations at an angle of attack of 10.5 degrees
            self.file10_sixth_AoA = "51.0_percent_semispan/10.5_degrees_AoA/51.0_percent_semispan_10_nodes_results_10.5_deg_AoA.csv"
            self.file17_sixth_AoA = "51.0_percent_semispan/10.5_degrees_AoA/51.0_percent_semispan_17_nodes_results_10.5_deg_AoA.csv"
            self.file25_sixth_AoA = "51.0_percent_semispan/10.5_degrees_AoA/51.0_percent_semispan_25_nodes_results_10.5_deg_AoA.csv"
            self.file45_sixth_AoA = "51.0_percent_semispan/10.5_degrees_AoA/51.0_percent_semispan_45_nodes_results_10.5_deg_AoA.csv"

            #Location of the experimental results
            file_experimental = "51.0_percent_semispan/51.0_percent_semispan_Experimental_Results_Wing_A.csv"

            self.x_LE = 24.99
            AoA = ['0', '2.1', '4.2', '6.3', '8.4', '10.5']

            labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
            title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.510$")    


        elif self.LE_Location == '65.3':

            #File locations at an angle of attack of 0 degrees
            self.file10_first_AoA = "65.3_percent_semispan/0_degrees_AoA/65.3_percent_semispan_10_nodes_results_0_deg_AoA.csv"
            self.file17_first_AoA = "65.3_percent_semispan/0_degrees_AoA/65.3_percent_semispan_17_nodes_results_0_deg_AoA.csv"
            self.file25_first_AoA = "65.3_percent_semispan/0_degrees_AoA/65.3_percent_semispan_25_nodes_results_0_deg_AoA.csv"
            self.file45_first_AoA = "65.3_percent_semispan/0_degrees_AoA/65.3_percent_semispan_45_nodes_results_0_deg_AoA.csv"

            #File locations at an angle of attack of 2.1 degrees
            self.file10_second_AoA = "65.3_percent_semispan/2.1_degrees_AoA/65.3_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
            self.file17_second_AoA = "65.3_percent_semispan/2.1_degrees_AoA/65.3_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
            self.file25_second_AoA = "65.3_percent_semispan/2.1_degrees_AoA/65.3_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
            self.file45_second_AoA = "65.3_percent_semispan/2.1_degrees_AoA/65.3_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

            #File locations at an angle of attack of 4.2 degrees
            self.file10_third_AoA = "65.3_percent_semispan/4.2_degrees_AoA/65.3_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
            self.file17_third_AoA = "65.3_percent_semispan/4.2_degrees_AoA/65.3_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
            self.file25_third_AoA = "65.3_percent_semispan/4.2_degrees_AoA/65.3_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
            self.file45_third_AoA = "65.3_percent_semispan/4.2_degrees_AoA/65.3_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

            #File locations at an angle of attack of 6.3 degrees
            self.file10_fourth_AoA = "65.3_percent_semispan/6.3_degrees_AoA/65.3_percent_semispan_10_nodes_results_6.3_deg_AoA.csv"
            self.file17_fourth_AoA = "65.3_percent_semispan/6.3_degrees_AoA/65.3_percent_semispan_17_nodes_results_6.3_deg_AoA.csv"
            self.file25_fourth_AoA = "65.3_percent_semispan/6.3_degrees_AoA/65.3_percent_semispan_25_nodes_results_6.3_deg_AoA.csv"
            self.file45_fourth_AoA = "65.3_percent_semispan/6.3_degrees_AoA/65.3_percent_semispan_45_nodes_results_6.3_deg_AoA.csv"

            #File locations at an angle of attack of 8.4 degrees
            self.file10_fifth_AoA = "65.3_percent_semispan/8.5_degrees_AoA/65.3_percent_semispan_10_nodes_results_8.5_deg_AoA.csv"
            self.file17_fifth_AoA = "65.3_percent_semispan/8.5_degrees_AoA/65.3_percent_semispan_17_nodes_results_8.5_deg_AoA.csv"
            self.file25_fifth_AoA = "65.3_percent_semispan/8.5_degrees_AoA/65.3_percent_semispan_25_nodes_results_8.5_deg_AoA.csv"
            self.file45_fifth_AoA = "65.3_percent_semispan/8.5_degrees_AoA/65.3_percent_semispan_45_nodes_results_8.5_deg_AoA.csv"

            #File locations at an angle of attack of 10.5 degrees
            self.file10_sixth_AoA = "65.3_percent_semispan/10.6_degrees_AoA/65.3_percent_semispan_10_nodes_results_10.6_deg_AoA.csv"
            self.file17_sixth_AoA = "65.3_percent_semispan/10.6_degrees_AoA/65.3_percent_semispan_17_nodes_results_10.6_deg_AoA.csv"
            self.file25_sixth_AoA = "65.3_percent_semispan/10.6_degrees_AoA/65.3_percent_semispan_25_nodes_results_10.6_deg_AoA.csv"
            self.file45_sixth_AoA = "65.3_percent_semispan/10.6_degrees_AoA/65.3_percent_semispan_45_nodes_results_10.6_deg_AoA.csv"

            #Location of the experimental results
            file_experimental = "51.0_percent_semispan/51.0_percent_semispan_Experimental_Results_Wing_A.csv"

            self.x_LE = 32
            AoA = ['0', '2.1', '4.2', '6.3', '8.5', '10.6']

            labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
            title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.653$")

        elif self.LE_Location == '89.8':
            #File locations at an angle of attack of 0 degrees
            #self.file10_first_AoA = "89.8_percent_semispan/0_degrees_AoA/89.8_percent_semispan_10_nodes_results_0_deg_AoA.csv"
            self.file17_first_AoA = "89.8_percent_semispan/0_degrees_AoA/89.8_percent_semispan_17_nodes_results_0_deg_AoA.csv"
            self.file25_first_AoA = "89.8_percent_semispan/0_degrees_AoA/89.8_percent_semispan_25_nodes_results_0_deg_AoA.csv"
            self.file45_first_AoA = "89.8_percent_semispan/0_degrees_AoA/89.8_percent_semispan_45_nodes_results_0_deg_AoA.csv"

            #File locations at an angle of attack of 4.3 degrees
            #self.file10_second_AoA = "89.8_percent_semispan/4.3_degrees_AoA/89.8_percent_semispan_10_nodes_results_4.3_deg_AoA.csv"
            self.file17_second_AoA = "89.8_percent_semispan/4.3_degrees_AoA/89.8_percent_semispan_17_nodes_results_4.3_deg_AoA.csv"
            self.file25_second_AoA = "89.8_percent_semispan/4.3_degrees_AoA/89.8_percent_semispan_25_nodes_results_4.3_deg_AoA.csv"
            self.file45_second_AoA = "89.8_percent_semispan/4.3_degrees_AoA/89.8_percent_semispan_45_nodes_results_4.3_deg_AoA.csv"

            #File locations at an angle of attack of 10.7 degrees
            #self.file10_third_AoA = "89.8_percent_semispan/10.7_degrees_AoA/89.8_percent_semispan_10_nodes_results_10.7_deg_AoA.csv"
            self.file17_third_AoA = "89.8_percent_semispan/10.7_degrees_AoA/89.8_percent_semispan_17_nodes_results_10.7_deg_AoA.csv"
            self.file25_third_AoA = "89.8_percent_semispan/10.7_degrees_AoA/89.8_percent_semispan_25_nodes_results_10.7_deg_AoA.csv"
            self.file45_third_AoA = "89.8_percent_semispan/10.7_degrees_AoA/89.8_percent_semispan_45_nodes_results_10.7_deg_AoA.csv"

            #Location of the experimental results
            file_experimental = "89.8_percent_semispan/89.8_percent_semispan_Experimental_Results_Wing_A.csv"
            
            self.x_LE = 44
            AoA = ['0', '4.3', '10.7']
            labels = ["17 Nodes", "25 Nodes", "45 Nodes"]
            #labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
            title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.898$")


        elif self.LE_Location == '94.9':
            #File locations at an angle of attack of 0 degrees
            #self.file10_first_AoA = "94.9_percent_semispan/0_degrees_AoA/94.9_percent_semispan_10_nodes_results_0_deg_AoA.csv"
            #self.file17_first_AoA = "94.9_percent_semispan/0_degrees_AoA/94.9_percent_semispan_17_nodes_results_0_deg_AoA.csv"
            #self.file25_first_AoA = "94.9_percent_semispan/0_degrees_AoA/94.9_percent_semispan_25_nodes_results_0_deg_AoA.csv"
            self.file45_first_AoA = "94.9_percent_semispan/0_degrees_AoA/94.9_percent_semispan_45_nodes_results_0_deg_AoA.csv"

            #File locations at an angle of attack of 2.1 degrees
            #self.file10_second_AoA = "94.9_percent_semispan/2.1_degrees_AoA/94.9_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
            #self.file17_second_AoA = "94.9_percent_semispan/2.1_degrees_AoA/94.9_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
            #self.file25_second_AoA = "94.9_percent_semispan/2.1_degrees_AoA/94.9_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
            self.file45_second_AoA = "94.9_percent_semispan/2.1_degrees_AoA/94.9_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

            #File locations at an angle of attack of 4.3 degrees
            #self.file10_third_AoA = "94.9_percent_semispan/4.3_degrees_AoA/94.9_percent_semispan_10_nodes_results_4.3_deg_AoA.csv"
            #self.file17_third_AoA = "94.9_percent_semispan/4.3_degrees_AoA/94.9_percent_semispan_17_nodes_results_4.3_deg_AoA.csv"
            #self.file25_third_AoA = "94.9_percent_semispan/4.3_degrees_AoA/94.9_percent_semispan_25_nodes_results_4.3_deg_AoA.csv"
            self.file45_third_AoA = "94.9_percent_semispan/4.3_degrees_AoA/94.9_percent_semispan_45_nodes_results_4.3_deg_AoA.csv"

            #File locations at an angle of attack of 6.4 degrees
            #self.file10_fourth_AoA = "94.9_percent_semispan/6.4_degrees_AoA/94.9_percent_semispan_10_nodes_results_6.4_deg_AoA.csv"
            #self.file17_fourth_AoA = "94.9_percent_semispan/6.4_degrees_AoA/94.9_percent_semispan_17_nodes_results_6.4_deg_AoA.csv"
            #self.file25_fourth_AoA = "94.9_percent_semispan/6.4_degrees_AoA/94.9_percent_semispan_25_nodes_results_6.4_deg_AoA.csv"
            self.file45_fourth_AoA = "94.9_percent_semispan/6.4_degrees_AoA/94.9_percent_semispan_45_nodes_results_6.4_deg_AoA.csv"

            #File locations at an angle of attack of 8.6 degrees
            #self.file10_fifth_AoA = "94.9_percent_semispan/8.6_degrees_AoA/94.9_percent_semispan_10_nodes_results_8.6_deg_AoA.csv"
            #self.file17_fifth_AoA = "94.9_percent_semispan/8.6_degrees_AoA/94.9_percent_semispan_17_nodes_results_8.6_deg_AoA.csv"
            #self.file25_fifth_AoA = "94.9_percent_semispan/8.6_degrees_AoA/94.9_percent_semispan_25_nodes_results_8.6_deg_AoA.csv"
            self.file45_fifth_AoA = "94.9_percent_semispan/8.6_degrees_AoA/94.9_percent_semispan_45_nodes_results_8.6_deg_AoA.csv"

            #File locations at an angle of attack of 10.7 degrees
            #self.file10_sixth_AoA = "94.9_percent_semispan/10.7_degrees_AoA/94.9_percent_semispan_10_nodes_results_10.7_deg_AoA.csv"
            #self.file17_sixth_AoA = "94.9_percent_semispan/10.7_degrees_AoA/94.9_percent_semispan_17_nodes_results_10.7_deg_AoA.csv"
            #self.file25_sixth_AoA = "94.9_percent_semispan/10.7_degrees_AoA/94.9_percent_semispan_25_nodes_results_10.7_deg_AoA.csv"
            self.file45_sixth_AoA = "94.9_percent_semispan/10.7_degrees_AoA/94.9_percent_semispan_45_nodes_results_10.7_deg_AoA.csv"

            #Location of the experimental results
            file_experimental = "94.9_percent_semispan/94.9_percent_semispan_Experimental_Results_Wing_A.csv"
            self.x_LE = 46.5

            labels =["45 Nodes"]
            AoA = ["0", "2.1", "4.3", "6.4", "8.6", "10.7"]
            #labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
            title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.949$")


        else:
            print("\n****************\nInvalid Entry. Please run script again.\n****************\n")
            quit()



        #Combine file locations into an iterable list
        if self.LE_Location == '89.8':
            file = [[self.file17_first_AoA, self.file25_first_AoA, self.file45_first_AoA],
                [self.file17_second_AoA, self.file25_second_AoA, self.file45_second_AoA],
                [self.file17_third_AoA, self.file25_third_AoA, self.file45_third_AoA]]

        elif self.LE_Location == '94.9':
            file = [[self.file45_first_AoA],
                [self.file45_second_AoA],
                [self.file45_third_AoA],
                [self.file45_fourth_AoA],
                [self.file45_fifth_AoA],
                [self.file45_sixth_AoA]]

        else:
            file = [[self.file10_first_AoA, self.file17_first_AoA, self.file25_first_AoA, self.file45_first_AoA],
                [self.file10_second_AoA, self.file17_second_AoA, self.file25_second_AoA, self.file45_second_AoA],
                [self.file10_third_AoA, self.file17_third_AoA, self.file25_third_AoA, self.file45_third_AoA],
                [self.file10_fourth_AoA, self.file17_fourth_AoA, self.file25_fourth_AoA, self.file45_fourth_AoA],
                [self.file10_fifth_AoA, self.file17_fifth_AoA, self.file25_fifth_AoA, self.file45_fifth_AoA],
                [self.file10_sixth_AoA, self.file17_sixth_AoA, self.file25_sixth_AoA, self.file45_sixth_AoA]]

        #Pull in experimental results
        Experimental = np.genfromtxt(file_experimental, delimiter=",", skip_header=2, dtype=float)


        #Identify Notes to include on plots
        y = r'$\alpha$ = '
        Notes = []
        for i in range(len(AoA)):
            Notes.append(y + AoA[i] + r'$^{\circ}$')

        #===Plotting Section===
        list_of_files = []
        #Iterate over all angles of attack
        for i in range(len(file)):

            complete_title = title + ', ' + Notes[i]
            #Iterate over all node densities
            for j in range(len(file[0])):

                self.Pressure_Plot(file[i][j], labels[j], self.LE_Location, chord)

            #Pull in experimental results for comparison at each angle of attack differentiating upper and lower surfaces
            upper_surface_count = Experimental[:,0].size//2 + Experimental[:,0].size%2

            plt.plot(Experimental[:upper_surface_count+1,0], Experimental[:upper_surface_count+1,i+1], ".",color="k", label="Exerimental Upper Surface")
            plt.plot(Experimental[upper_surface_count:,0], Experimental[upper_surface_count:,i+1], "*",color="k", label="Exerimental Lower Surface")



            #Plot the figure containing all curves
            plt.title(complete_title)
            plt.xlabel('x/c')
            plt.ylabel(r"$C_p$")
            plt.gca().invert_yaxis()
            if self.LE_Location == '94.9':
                plt.legend()
            else:
                plt.legend(ncol = 2)

            #Save the figure in its appropriate location
            if json_vals["plots"]["save plots"]:
                filename = self.LE_Location + "_percent_semispan/plots_" + self.LE_Location + "_percent_semispan/" + json_vals["plots"]["save plot type"] + "_plots/" + AoA[i] + "degrees_AoA_plot." + json_vals["plots"]["save plot type"]
                list_of_files.append(filename)
                plt.savefig(filename)

            plt.show()


        # Combine all pdf output files for each case in respective folder locations
        # merger = PdfFileMerger()
        # Combined_name = self.LE_Location + '_percent_semispan/plots_' + self.LE_Location + '_percent_semispan/Combined_plots_' + self.LE_Location + '_percent_semispan.pdf'
        # Combined_name_copy = 'Plot_Summary/' + 'Combined_plots_' + self.LE_Location + '_percent_semispan.pdf'

        # for k in range(len(list_of_files)):
        #     merger.append(list_of_files[k])

        # merger.write(Combined_name)
        # merger.write(Combined_name_copy)
        # merger.close()

inputfile = "Swept_half_wing_conditions_input.json"

json_string = open(inputfile).read()
json_vals = json.loads(json_string)

if json_vals["plots"]["show_all"]:
    for i in range(len(json_vals["geometry"]["Percent Semispan Locations"])):
        Swept_Plotting(json_vals["geometry"]["Percent Semispan Locations"][i]).get_data()
else:
    x = input("Enter percent semispan for analysis results. Options are 4.1, 8.2, 16.3, 24.5, 36.7, 51.0, 65.3, 89.8, 94.9:   ",)
    Swept_Plotting(x).get_data()
