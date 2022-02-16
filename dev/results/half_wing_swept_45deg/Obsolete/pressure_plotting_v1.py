#This file is to take .csv or .txt files of lift distributions and plot
#them allowing comparison to experimental results

import numpy as np
import matplotlib.pyplot as plt
#from PyPDF2 import PdfFileMerger

def Pressure_Plot(file, Plot_label, LE_loc):
    Data = np.genfromtxt(file, delimiter=",", skip_header=1, dtype=float)

    C_p = Data[:,0]
    x = Data[:,1]

    # Ensure the plots are colored consistently
    if LE_loc == '89.8':
            if file == file17_first_AoA or file == file17_second_AoA or file == file17_third_AoA:

                colors = 'b'

            elif file == file25_first_AoA or file == file25_second_AoA or file == file25_third_AoA:

                colors = 'c'

            elif file == file45_first_AoA or file == file45_second_AoA or file == file45_third_AoA:

                colors = 'r'

            else:
                colors = 'y'
            
    elif LE_loc == '94.9':
            if file == file45_first_AoA or file == file45_second_AoA or file == file45_third_AoA or file == file45_fourth_AoA or file == file45_fifth_AoA or file == file45_sixth_AoA:

                colors = 'r'

            else:
                colors = 'y'

    else:
        if file == file10_first_AoA or file == file10_second_AoA or file == file10_third_AoA or file == file10_fourth_AoA or file == file10_fifth_AoA or file == file10_sixth_AoA:

            colors = 'g'

        elif file == file17_first_AoA or file == file17_second_AoA or file == file17_third_AoA or file == file17_fourth_AoA or file == file17_fifth_AoA or file == file17_sixth_AoA:

            colors = 'b'

        elif file == file25_first_AoA or file == file25_second_AoA or file == file25_third_AoA or file == file25_fourth_AoA or file == file25_fifth_AoA or file == file25_sixth_AoA:

            colors = 'c'

        elif file == file45_first_AoA or file == file45_second_AoA or file == file45_third_AoA or file == file45_fourth_AoA or file == file45_fifth_AoA or file == file45_sixth_AoA:

            colors = 'r'

        else:
            colors = 'y'
    
    plt.plot((x-x_LE)/chord,C_p, label=Plot_label, color = colors)

    return

#Main

LE_Location = input("Enter percent semispan for analysis results. Options are 4.1, 8.2, 16.3, 24.5, 36.7, 51.0, 65.3, 89.8, 94.9:   ",)
chord = 20

#Pull In Results Depending on Location Input
#In addition to defining leading edge location and respective labels for plot
if LE_Location == '4.1':
    #File locations at an angle of attack of 0 degrees   
    file10_first_AoA = "4.1_percent_semispan/0degrees_AoA/4.1_percent_semispan_10_nodes_results.csv"
    file17_first_AoA = "4.1_percent_semispan/0degrees_AoA/4.1_percent_semispan_17_nodes_results.csv"
    file25_first_AoA = "4.1_percent_semispan/0degrees_AoA/4.1_percent_semispan_25_nodes_results.csv"
    file45_first_AoA = "4.1_percent_semispan/0degrees_AoA/4.1_percent_semispan_45_nodes_results.csv"
    
    #File locations at an angle of attack of 2.1 degrees
    file10_second_AoA = "4.1_percent_semispan/2.1_degrees_AoA/4.1_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
    file17_second_AoA = "4.1_percent_semispan/2.1_degrees_AoA/4.1_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
    file25_second_AoA = "4.1_percent_semispan/2.1_degrees_AoA/4.1_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
    file45_second_AoA = "4.1_percent_semispan/2.1_degrees_AoA/4.1_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

    #File locations at an angle of attack of 4.2 degrees
    file10_third_AoA = "4.1_percent_semispan/4.2_degrees_AoA/4.1_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
    file17_third_AoA = "4.1_percent_semispan/4.2_degrees_AoA/4.1_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
    file25_third_AoA = "4.1_percent_semispan/4.2_degrees_AoA/4.1_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
    file45_third_AoA = "4.1_percent_semispan/4.2_degrees_AoA/4.1_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

    #File locations at an angle of attack of 6.2 degrees
    file10_fourth_AoA = "4.1_percent_semispan/6.2_degrees_AoA/4.1_percent_semispan_10_nodes_results_6.2_deg_AoA.csv"
    file17_fourth_AoA = "4.1_percent_semispan/6.2_degrees_AoA/4.1_percent_semispan_17_nodes_results_6.2_deg_AoA.csv"
    file25_fourth_AoA = "4.1_percent_semispan/6.2_degrees_AoA/4.1_percent_semispan_25_nodes_results_6.2_deg_AoA.csv"
    file45_fourth_AoA = "4.1_percent_semispan/6.2_degrees_AoA/4.1_percent_semispan_45_nodes_results_6.2_deg_AoA.csv"

    #File locations at an angle of attack of 8.3 degrees
    file10_fifth_AoA = "4.1_percent_semispan/8.3_degrees_AoA/4.1_percent_semispan_10_nodes_results_8.3_deg_AoA.csv"
    file17_fifth_AoA = "4.1_percent_semispan/8.3_degrees_AoA/4.1_percent_semispan_17_nodes_results_8.3_deg_AoA.csv"
    file25_fifth_AoA = "4.1_percent_semispan/8.3_degrees_AoA/4.1_percent_semispan_25_nodes_results_8.3_deg_AoA.csv"
    file45_fifth_AoA = "4.1_percent_semispan/8.3_degrees_AoA/4.1_percent_semispan_45_nodes_results_8.3_deg_AoA.csv"

    #File locations at an angle of attack of 10.4 degrees
    file10_sixth_AoA = "4.1_percent_semispan/10.4_degrees_AoA/4.1_percent_semispan_10_nodes_results_10.4_deg_AoA.csv"
    file17_sixth_AoA = "4.1_percent_semispan/10.4_degrees_AoA/4.1_percent_semispan_17_nodes_results_10.4_deg_AoA.csv"
    file25_sixth_AoA = "4.1_percent_semispan/10.4_degrees_AoA/4.1_percent_semispan_25_nodes_results_10.4_deg_AoA.csv"
    file45_sixth_AoA = "4.1_percent_semispan/10.4_degrees_AoA/4.1_percent_semispan_45_nodes_results_10.4_deg_AoA.csv"

    #Location of the experimental results
    file_experimental = "4.1_percent_semispan/4.1_percent_semispan_Experimental_Results_Wing_A.csv"
    
    x_LE = 2
    AoA = ['0', '2.1', '4.2', '6.2', '8.3', '10.4']
    
    #Identify labels and notes to include on plots
    labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
    title = r"Pressure Coefficient at $\frac{2y}{b} = 0.041$"
        

elif LE_Location == '8.2':
    #File locations at an angle of attack of 0 degrees
    file10_first_AoA = "8.2_percent_semispan/0degrees_AoA/8.2_percent_semispan_10_nodes_results.csv"
    file17_first_AoA = "8.2_percent_semispan/0degrees_AoA/8.2_percent_semispan_17_nodes_results.csv"
    file25_first_AoA = "8.2_percent_semispan/0degrees_AoA/8.2_percent_semispan_25_nodes_results.csv"
    file45_first_AoA = "8.2_percent_semispan/0degrees_AoA/8.2_percent_semispan_45_nodes_results.csv"

    #File locations at an angle of attack of 2.1 degrees
    file10_second_AoA = "8.2_percent_semispan/2.1_degrees_AoA/8.2_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
    file17_second_AoA = "8.2_percent_semispan/2.1_degrees_AoA/8.2_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
    file25_second_AoA = "8.2_percent_semispan/2.1_degrees_AoA/8.2_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
    file45_second_AoA = "8.2_percent_semispan/2.1_degrees_AoA/8.2_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

    #File locations at an angle of attack of 4.2 degrees
    file10_third_AoA = "8.2_percent_semispan/4.2_degrees_AoA/8.2_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
    file17_third_AoA = "8.2_percent_semispan/4.2_degrees_AoA/8.2_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
    file25_third_AoA = "8.2_percent_semispan/4.2_degrees_AoA/8.2_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
    file45_third_AoA = "8.2_percent_semispan/4.2_degrees_AoA/8.2_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

    #File locations at an angle of attack of 6.2 degrees
    file10_fourth_AoA = "8.2_percent_semispan/6.2_degrees_AoA/8.2_percent_semispan_10_nodes_results_6.2_deg_AoA.csv"
    file17_fourth_AoA = "8.2_percent_semispan/6.2_degrees_AoA/8.2_percent_semispan_17_nodes_results_6.2_deg_AoA.csv"
    file25_fourth_AoA = "8.2_percent_semispan/6.2_degrees_AoA/8.2_percent_semispan_25_nodes_results_6.2_deg_AoA.csv"
    file45_fourth_AoA = "8.2_percent_semispan/6.2_degrees_AoA/8.2_percent_semispan_45_nodes_results_6.2_deg_AoA.csv"

    #File locations at an angle of attack of 8.3 degrees
    file10_fifth_AoA = "8.2_percent_semispan/8.3_degrees_AoA/8.2_percent_semispan_10_nodes_results_8.3_deg_AoA.csv"
    file17_fifth_AoA = "8.2_percent_semispan/8.3_degrees_AoA/8.2_percent_semispan_17_nodes_results_8.3_deg_AoA.csv"
    file25_fifth_AoA = "8.2_percent_semispan/8.3_degrees_AoA/8.2_percent_semispan_25_nodes_results_8.3_deg_AoA.csv"
    file45_fifth_AoA = "8.2_percent_semispan/8.3_degrees_AoA/8.2_percent_semispan_45_nodes_results_8.3_deg_AoA.csv"

    #File locations at an angle of attack of 10.4 degrees
    file10_sixth_AoA = "8.2_percent_semispan/10.4_degrees_AoA/8.2_percent_semispan_10_nodes_results_10.4_deg_AoA.csv"
    file17_sixth_AoA = "8.2_percent_semispan/10.4_degrees_AoA/8.2_percent_semispan_17_nodes_results_10.4_deg_AoA.csv"
    file25_sixth_AoA = "8.2_percent_semispan/10.4_degrees_AoA/8.2_percent_semispan_25_nodes_results_10.4_deg_AoA.csv"
    file45_sixth_AoA = "8.2_percent_semispan/10.4_degrees_AoA/8.2_percent_semispan_45_nodes_results_10.4_deg_AoA.csv"

    #Location of the experimental results
    file_experimental = "8.2_percent_semispan/8.2_percent_semispan_Experimental_Results_Wing_A.csv"
   
    x_LE = 4
    AoA = ['0', '2.1', '4.2', '6.2', '8.3', '10.4']
    
    #Identify labels and title to include on plots
    labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
    title = r"Pressure Coefficient at $\frac{2y}{b} = 0.082$"


elif LE_Location == '16.3':
    #File locations at an angle of attack of 0 degrees
    file10_first_AoA = "16.3_percent_semispan/0degrees_AoA/16.3_percent_semispan_10_nodes_results.csv"
    file17_first_AoA = "16.3_percent_semispan/0degrees_AoA/16.3_percent_semispan_17_nodes_results.csv"
    file25_first_AoA = "16.3_percent_semispan/0degrees_AoA/16.3_percent_semispan_25_nodes_results.csv"
    file45_first_AoA = "16.3_percent_semispan/0degrees_AoA/16.3_percent_semispan_45_nodes_results.csv"
    
    #File locations at an angle of attack of 2.1 degrees
    file10_second_AoA = "16.3_percent_semispan/2.1_degrees_AoA/16.3_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
    file17_second_AoA = "16.3_percent_semispan/2.1_degrees_AoA/16.3_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
    file25_second_AoA = "16.3_percent_semispan/2.1_degrees_AoA/16.3_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
    file45_second_AoA = "16.3_percent_semispan/2.1_degrees_AoA/16.3_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

    #File locations at an angle of attack of 4.2 degrees
    file10_third_AoA = "16.3_percent_semispan/4.2_degrees_AoA/16.3_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
    file17_third_AoA = "16.3_percent_semispan/4.2_degrees_AoA/16.3_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
    file25_third_AoA = "16.3_percent_semispan/4.2_degrees_AoA/16.3_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
    file45_third_AoA = "16.3_percent_semispan/4.2_degrees_AoA/16.3_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

    #File locations at an angle of attack of 6.2 degrees
    file10_fourth_AoA = "16.3_percent_semispan/6.2_degrees_AoA/16.3_percent_semispan_10_nodes_results_6.2_deg_AoA.csv"
    file17_fourth_AoA = "16.3_percent_semispan/6.2_degrees_AoA/16.3_percent_semispan_17_nodes_results_6.2_deg_AoA.csv"
    file25_fourth_AoA = "16.3_percent_semispan/6.2_degrees_AoA/16.3_percent_semispan_25_nodes_results_6.2_deg_AoA.csv"
    file45_fourth_AoA = "16.3_percent_semispan/6.2_degrees_AoA/16.3_percent_semispan_45_nodes_results_6.2_deg_AoA.csv"

    #File locations at an angle of attack of 8.3 degrees
    file10_fifth_AoA = "16.3_percent_semispan/8.3_degrees_AoA/16.3_percent_semispan_10_nodes_results_8.3_deg_AoA.csv"
    file17_fifth_AoA = "16.3_percent_semispan/8.3_degrees_AoA/16.3_percent_semispan_17_nodes_results_8.3_deg_AoA.csv"
    file25_fifth_AoA = "16.3_percent_semispan/8.3_degrees_AoA/16.3_percent_semispan_25_nodes_results_8.3_deg_AoA.csv"
    file45_fifth_AoA = "16.3_percent_semispan/8.3_degrees_AoA/16.3_percent_semispan_45_nodes_results_8.3_deg_AoA.csv"

    #File locations at an angle of attack of 10.4 degrees
    file10_sixth_AoA = "16.3_percent_semispan/10.4_degrees_AoA/16.3_percent_semispan_10_nodes_results_10.4_deg_AoA.csv"
    file17_sixth_AoA = "16.3_percent_semispan/10.4_degrees_AoA/16.3_percent_semispan_17_nodes_results_10.4_deg_AoA.csv"
    file25_sixth_AoA = "16.3_percent_semispan/10.4_degrees_AoA/16.3_percent_semispan_25_nodes_results_10.4_deg_AoA.csv"
    file45_sixth_AoA = "16.3_percent_semispan/10.4_degrees_AoA/16.3_percent_semispan_45_nodes_results_10.4_deg_AoA.csv"

    #Location of the experimental results
    file_experimental = "16.3_percent_semispan/16.3_percent_semispan_Experimental_Results_Wing_A.csv"
    
    x_LE = 7.99
    AoA = ['0', '2.1', '4.2', '6.2', '8.3', '10.4']
    
    #Identify labels and title to include on plots
    labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
    title = r"Pressure Coefficient at $\frac{2y}{b} = 0.163$"


elif LE_Location == '24.5':
    #File locations at an angle of attack of 0 degrees
    file10_first_AoA = "24.5_percent_semispan/0degrees_AoA/24.5_percent_semispan_10_nodes_results.csv"
    file17_first_AoA = "24.5_percent_semispan/0degrees_AoA/24.5_percent_semispan_17_nodes_results.csv"
    file25_first_AoA = "24.5_percent_semispan/0degrees_AoA/24.5_percent_semispan_25_nodes_results.csv"
    file45_first_AoA = "24.5_percent_semispan/0degrees_AoA/24.5_percent_semispan_45_nodes_results.csv"

    #File locations at an angle of attack of 2.1 degrees
    file10_second_AoA = "24.5_percent_semispan/2.1_degrees_AoA/24.5_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
    file17_second_AoA = "24.5_percent_semispan/2.1_degrees_AoA/24.5_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
    file25_second_AoA = "24.5_percent_semispan/2.1_degrees_AoA/24.5_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
    file45_second_AoA = "24.5_percent_semispan/2.1_degrees_AoA/24.5_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

    #File locations at an angle of attack of 4.2 degrees
    file10_third_AoA = "24.5_percent_semispan/4.2_degrees_AoA/24.5_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
    file17_third_AoA = "24.5_percent_semispan/4.2_degrees_AoA/24.5_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
    file25_third_AoA = "24.5_percent_semispan/4.2_degrees_AoA/24.5_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
    file45_third_AoA = "24.5_percent_semispan/4.2_degrees_AoA/24.5_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

    #File locations at an angle of attack of 6.2 degrees
    file10_fourth_AoA = "24.5_percent_semispan/6.2_degrees_AoA/24.5_percent_semispan_10_nodes_results_6.2_deg_AoA.csv"
    file17_fourth_AoA = "24.5_percent_semispan/6.2_degrees_AoA/24.5_percent_semispan_17_nodes_results_6.2_deg_AoA.csv"
    file25_fourth_AoA = "24.5_percent_semispan/6.2_degrees_AoA/24.5_percent_semispan_25_nodes_results_6.2_deg_AoA.csv"
    file45_fourth_AoA = "24.5_percent_semispan/6.2_degrees_AoA/24.5_percent_semispan_45_nodes_results_6.2_deg_AoA.csv"

    #File locations at an angle of attack of 8.3 degrees
    file10_fifth_AoA = "24.5_percent_semispan/8.3_degrees_AoA/24.5_percent_semispan_10_nodes_results_8.3_deg_AoA.csv"
    file17_fifth_AoA = "24.5_percent_semispan/8.3_degrees_AoA/24.5_percent_semispan_17_nodes_results_8.3_deg_AoA.csv"
    file25_fifth_AoA = "24.5_percent_semispan/8.3_degrees_AoA/24.5_percent_semispan_25_nodes_results_8.3_deg_AoA.csv"
    file45_fifth_AoA = "24.5_percent_semispan/8.3_degrees_AoA/24.5_percent_semispan_45_nodes_results_8.3_deg_AoA.csv"

    #File locations at an angle of attack of 10.4 degrees
    file10_sixth_AoA = "24.5_percent_semispan/10.4_degrees_AoA/24.5_percent_semispan_10_nodes_results_10.4_deg_AoA.csv"
    file17_sixth_AoA = "24.5_percent_semispan/10.4_degrees_AoA/24.5_percent_semispan_17_nodes_results_10.4_deg_AoA.csv"
    file25_sixth_AoA = "24.5_percent_semispan/10.4_degrees_AoA/24.5_percent_semispan_25_nodes_results_10.4_deg_AoA.csv"
    file45_sixth_AoA = "24.5_percent_semispan/10.4_degrees_AoA/24.5_percent_semispan_45_nodes_results_10.4_deg_AoA.csv"

    #Location of the experimental results
    file_experimental = "24.5_percent_semispan/24.5_percent_semispan_Experimental_Results_Wing_A.csv"

    x_LE = 12
    AoA = ['0', '2.1', '4.2', '6.2', '8.3', '10.4']

    labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
    title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.245$")


elif LE_Location == '36.7':
    #File locations at an angle of attack of 0 degrees
    file10_first_AoA = "36.7_percent_semispan/0degrees_AoA/36.7_percent_semispan_10_nodes_results.csv"
    file17_first_AoA = "36.7_percent_semispan/0degrees_AoA/36.7_percent_semispan_17_nodes_results.csv"
    file25_first_AoA = "36.7_percent_semispan/0degrees_AoA/36.7_percent_semispan_25_nodes_results.csv"
    file45_first_AoA = "36.7_percent_semispan/0degrees_AoA/36.7_percent_semispan_45_nodes_results.csv"

    #File locations at an angle of attack of 2.1 degrees
    file10_second_AoA = "36.7_percent_semispan/2.1_degrees_AoA/36.7_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
    file17_second_AoA = "36.7_percent_semispan/2.1_degrees_AoA/36.7_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
    file25_second_AoA = "36.7_percent_semispan/2.1_degrees_AoA/36.7_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
    file45_second_AoA = "36.7_percent_semispan/2.1_degrees_AoA/36.7_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

    #File locations at an angle of attack of 4.2 degrees
    file10_third_AoA = "36.7_percent_semispan/4.2_degrees_AoA/36.7_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
    file17_third_AoA = "36.7_percent_semispan/4.2_degrees_AoA/36.7_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
    file25_third_AoA = "36.7_percent_semispan/4.2_degrees_AoA/36.7_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
    file45_third_AoA = "36.7_percent_semispan/4.2_degrees_AoA/36.7_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

    #File locations at an angle of attack of 6.3 degrees
    file10_fourth_AoA = "36.7_percent_semispan/6.3_degrees_AoA/36.7_percent_semispan_10_nodes_results_6.3_deg_AoA.csv"
    file17_fourth_AoA = "36.7_percent_semispan/6.3_degrees_AoA/36.7_percent_semispan_17_nodes_results_6.3_deg_AoA.csv"
    file25_fourth_AoA = "36.7_percent_semispan/6.3_degrees_AoA/36.7_percent_semispan_25_nodes_results_6.3_deg_AoA.csv"
    file45_fourth_AoA = "36.7_percent_semispan/6.3_degrees_AoA/36.7_percent_semispan_45_nodes_results_6.3_deg_AoA.csv"

    #File locations at an angle of attack of 8.4 degrees
    file10_fifth_AoA = "36.7_percent_semispan/8.4_degrees_AoA/36.7_percent_semispan_10_nodes_results_8.4_deg_AoA.csv"
    file17_fifth_AoA = "36.7_percent_semispan/8.4_degrees_AoA/36.7_percent_semispan_17_nodes_results_8.4_deg_AoA.csv"
    file25_fifth_AoA = "36.7_percent_semispan/8.4_degrees_AoA/36.7_percent_semispan_25_nodes_results_8.4_deg_AoA.csv"
    file45_fifth_AoA = "36.7_percent_semispan/8.4_degrees_AoA/36.7_percent_semispan_45_nodes_results_8.4_deg_AoA.csv"

    #File locations at an angle of attack of 10.5 degrees
    file10_sixth_AoA = "36.7_percent_semispan/10.5_degrees_AoA/36.7_percent_semispan_10_nodes_results_10.5_deg_AoA.csv"
    file17_sixth_AoA = "36.7_percent_semispan/10.5_degrees_AoA/36.7_percent_semispan_17_nodes_results_10.5_deg_AoA.csv"
    file25_sixth_AoA = "36.7_percent_semispan/10.5_degrees_AoA/36.7_percent_semispan_25_nodes_results_10.5_deg_AoA.csv"
    file45_sixth_AoA = "36.7_percent_semispan/10.5_degrees_AoA/36.7_percent_semispan_45_nodes_results_10.5_deg_AoA.csv"

    #Location of the experimental results
    file_experimental = "36.7_percent_semispan/36.7_percent_semispan_Experimental_Results_Wing_A.csv"

    x_LE = 17.98
    AoA = ['0', '2.1', '4.2', '6.3', '8.4', '10.5']

    labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
    title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.367$")


elif LE_Location == '51.0' or LE_Location == '51':
    LE_Location = "51.0"

    #File locations at an angle of attack of 0 degrees
    file10_first_AoA = "51.0_percent_semispan/0degrees_AoA/51.0_percent_semispan_10_nodes_results.csv"
    file17_first_AoA = "51.0_percent_semispan/0degrees_AoA/51.0_percent_semispan_17_nodes_results.csv"
    file25_first_AoA = "51.0_percent_semispan/0degrees_AoA/51.0_percent_semispan_25_nodes_results.csv"
    file45_first_AoA = "51.0_percent_semispan/0degrees_AoA/51.0_percent_semispan_45_nodes_results.csv"

    #File locations at an angle of attack of 2.1 degrees
    file10_second_AoA = "51.0_percent_semispan/2.1_degrees_AoA/51.0_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
    file17_second_AoA = "51.0_percent_semispan/2.1_degrees_AoA/51.0_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
    file25_second_AoA = "51.0_percent_semispan/2.1_degrees_AoA/51.0_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
    file45_second_AoA = "51.0_percent_semispan/2.1_degrees_AoA/51.0_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

    #File locations at an angle of attack of 4.2 degrees
    file10_third_AoA = "51.0_percent_semispan/4.2_degrees_AoA/51.0_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
    file17_third_AoA = "51.0_percent_semispan/4.2_degrees_AoA/51.0_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
    file25_third_AoA = "51.0_percent_semispan/4.2_degrees_AoA/51.0_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
    file45_third_AoA = "51.0_percent_semispan/4.2_degrees_AoA/51.0_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

    #File locations at an angle of attack of 6.3 degrees
    file10_fourth_AoA = "51.0_percent_semispan/6.3_degrees_AoA/51.0_percent_semispan_10_nodes_results_6.3_deg_AoA.csv"
    file17_fourth_AoA = "51.0_percent_semispan/6.3_degrees_AoA/51.0_percent_semispan_17_nodes_results_6.3_deg_AoA.csv"
    file25_fourth_AoA = "51.0_percent_semispan/6.3_degrees_AoA/51.0_percent_semispan_25_nodes_results_6.3_deg_AoA.csv"
    file45_fourth_AoA = "51.0_percent_semispan/6.3_degrees_AoA/51.0_percent_semispan_45_nodes_results_6.3_deg_AoA.csv"

    #File locations at an angle of attack of 8.4 degrees
    file10_fifth_AoA = "51.0_percent_semispan/8.4_degrees_AoA/51.0_percent_semispan_10_nodes_results_8.4_deg_AoA.csv"
    file17_fifth_AoA = "51.0_percent_semispan/8.4_degrees_AoA/51.0_percent_semispan_17_nodes_results_8.4_deg_AoA.csv"
    file25_fifth_AoA = "51.0_percent_semispan/8.4_degrees_AoA/51.0_percent_semispan_25_nodes_results_8.4_deg_AoA.csv"
    file45_fifth_AoA = "51.0_percent_semispan/8.4_degrees_AoA/51.0_percent_semispan_45_nodes_results_8.4_deg_AoA.csv"

    #File locations at an angle of attack of 10.5 degrees
    file10_sixth_AoA = "51.0_percent_semispan/10.5_degrees_AoA/51.0_percent_semispan_10_nodes_results_10.5_deg_AoA.csv"
    file17_sixth_AoA = "51.0_percent_semispan/10.5_degrees_AoA/51.0_percent_semispan_17_nodes_results_10.5_deg_AoA.csv"
    file25_sixth_AoA = "51.0_percent_semispan/10.5_degrees_AoA/51.0_percent_semispan_25_nodes_results_10.5_deg_AoA.csv"
    file45_sixth_AoA = "51.0_percent_semispan/10.5_degrees_AoA/51.0_percent_semispan_45_nodes_results_10.5_deg_AoA.csv"

    #Location of the experimental results
    file_experimental = "51.0_percent_semispan/51.0_percent_semispan_Experimental_Results_Wing_A.csv"

    x_LE = 24.99
    AoA = ['0', '2.1', '4.2', '6.3', '8.4', '10.5']

    labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
    title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.510$")    


elif LE_Location == '65.3':
    #File locations at an angle of attack of 0 degrees
    file10_first_AoA = "65.3_percent_semispan/0degrees_AoA/65.3_percent_semispan_10_nodes_results.csv"
    file17_first_AoA = "65.3_percent_semispan/0degrees_AoA/65.3_percent_semispan_17_nodes_results.csv"
    file25_first_AoA = "65.3_percent_semispan/0degrees_AoA/65.3_percent_semispan_25_nodes_results.csv"
    file45_first_AoA = "65.3_percent_semispan/0degrees_AoA/65.3_percent_semispan_45_nodes_results.csv"

    #File locations at an angle of attack of 2.1 degrees
    file10_second_AoA = "65.3_percent_semispan/2.1_degrees_AoA/65.3_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
    file17_second_AoA = "65.3_percent_semispan/2.1_degrees_AoA/65.3_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
    file25_second_AoA = "65.3_percent_semispan/2.1_degrees_AoA/65.3_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
    file45_second_AoA = "65.3_percent_semispan/2.1_degrees_AoA/65.3_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

    #File locations at an angle of attack of 4.2 degrees
    file10_third_AoA = "65.3_percent_semispan/4.2_degrees_AoA/65.3_percent_semispan_10_nodes_results_4.2_deg_AoA.csv"
    file17_third_AoA = "65.3_percent_semispan/4.2_degrees_AoA/65.3_percent_semispan_17_nodes_results_4.2_deg_AoA.csv"
    file25_third_AoA = "65.3_percent_semispan/4.2_degrees_AoA/65.3_percent_semispan_25_nodes_results_4.2_deg_AoA.csv"
    file45_third_AoA = "65.3_percent_semispan/4.2_degrees_AoA/65.3_percent_semispan_45_nodes_results_4.2_deg_AoA.csv"

    #File locations at an angle of attack of 6.3 degrees
    file10_fourth_AoA = "65.3_percent_semispan/6.3_degrees_AoA/65.3_percent_semispan_10_nodes_results_6.3_deg_AoA.csv"
    file17_fourth_AoA = "65.3_percent_semispan/6.3_degrees_AoA/65.3_percent_semispan_17_nodes_results_6.3_deg_AoA.csv"
    file25_fourth_AoA = "65.3_percent_semispan/6.3_degrees_AoA/65.3_percent_semispan_25_nodes_results_6.3_deg_AoA.csv"
    file45_fourth_AoA = "65.3_percent_semispan/6.3_degrees_AoA/65.3_percent_semispan_45_nodes_results_6.3_deg_AoA.csv"

    #File locations at an angle of attack of 8.5 degrees
    file10_fifth_AoA = "65.3_percent_semispan/8.5_degrees_AoA/65.3_percent_semispan_10_nodes_results_8.5_deg_AoA.csv"
    file17_fifth_AoA = "65.3_percent_semispan/8.5_degrees_AoA/65.3_percent_semispan_17_nodes_results_8.5_deg_AoA.csv"
    file25_fifth_AoA = "65.3_percent_semispan/8.5_degrees_AoA/65.3_percent_semispan_25_nodes_results_8.5_deg_AoA.csv"
    file45_fifth_AoA = "65.3_percent_semispan/8.5_degrees_AoA/65.3_percent_semispan_45_nodes_results_8.5_deg_AoA.csv"

    #File locations at an angle of attack of 10.6 degrees
    file10_sixth_AoA = "65.3_percent_semispan/10.6_degrees_AoA/65.3_percent_semispan_10_nodes_results_10.6_deg_AoA.csv"
    file17_sixth_AoA = "65.3_percent_semispan/10.6_degrees_AoA/65.3_percent_semispan_17_nodes_results_10.6_deg_AoA.csv"
    file25_sixth_AoA = "65.3_percent_semispan/10.6_degrees_AoA/65.3_percent_semispan_25_nodes_results_10.6_deg_AoA.csv"
    file45_sixth_AoA = "65.3_percent_semispan/10.6_degrees_AoA/65.3_percent_semispan_45_nodes_results_10.6_deg_AoA.csv"

    #Location of the experimental results
    file_experimental = "65.3_percent_semispan/65.3_percent_semispan_Experimental_Results_Wing_A.csv"

    x_LE = 32
    AoA = ['0', '2.1', '4.2', '6.3', '8.5', '10.6']

    labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
    title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.653$")


elif LE_Location == '89.8':
    #File locations at an angle of attack of 0 degrees
    #file10_first_AoA = "89.8_percent_semispan/0degrees_AoA/89.8_percent_semispan_10_nodes_results.csv"
    file17_first_AoA = "89.8_percent_semispan/0degrees_AoA/89.8_percent_semispan_17_nodes_results.csv"
    file25_first_AoA = "89.8_percent_semispan/0degrees_AoA/89.8_percent_semispan_25_nodes_results.csv"
    file45_first_AoA = "89.8_percent_semispan/0degrees_AoA/89.8_percent_semispan_45_nodes_results.csv"

    #File locations at an angle of attack of 4.3 degrees
    #file10_second_AoA = "89.8_percent_semispan/4.3_degrees_AoA/89.8_percent_semispan_10_nodes_results_4.3_deg_AoA.csv"
    file17_second_AoA = "89.8_percent_semispan/4.3_degrees_AoA/89.8_percent_semispan_17_nodes_results_4.3_deg_AoA.csv"
    file25_second_AoA = "89.8_percent_semispan/4.3_degrees_AoA/89.8_percent_semispan_25_nodes_results_4.3_deg_AoA.csv"
    file45_second_AoA = "89.8_percent_semispan/4.3_degrees_AoA/89.8_percent_semispan_45_nodes_results_4.3_deg_AoA.csv"

    #File locations at an angle of attack of 10.7 degrees
    #file10_third_AoA = "89.8_percent_semispan/10.7_degrees_AoA/89.8_percent_semispan_10_nodes_results_10.7_deg_AoA.csv"
    file17_third_AoA = "89.8_percent_semispan/10.7_degrees_AoA/89.8_percent_semispan_17_nodes_results_10.7_deg_AoA.csv"
    file25_third_AoA = "89.8_percent_semispan/10.7_degrees_AoA/89.8_percent_semispan_25_nodes_results_10.7_deg_AoA.csv"
    file45_third_AoA = "89.8_percent_semispan/10.7_degrees_AoA/89.8_percent_semispan_45_nodes_results_10.7_deg_AoA.csv"

    #Location of the experimental results
    file_experimental = "89.8_percent_semispan/89.8_percent_semispan_Experimental_Results_Wing_A.csv"
    
    x_LE = 44
    AoA = ['0', '4.3', '10.7']
    labels = ["17 Nodes", "25 Nodes", "45 Nodes"]
    #labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
    title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.898$")


elif LE_Location == '94.9':
    #File locations at an angle of attack of 0 degrees
    #file10_first_AoA = "94.9_percent_semispan/0degrees_AoA/94.9_percent_semispan_10_nodes_results.csv"
    #file17_first_AoA = "94.9_percent_semispan/0degrees_AoA/94.9_percent_semispan_17_nodes_results.csv"
    #file25_first_AoA = "94.9_percent_semispan/0degrees_AoA/94.9_percent_semispan_25_nodes_results.csv"
    file45_first_AoA = "94.9_percent_semispan/0degrees_AoA/94.9_percent_semispan_45_nodes_results.csv"

    #File locations at an angle of attack of 2.1 degrees
    #file10_second_AoA = "94.9_percent_semispan/2.1_degrees_AoA/94.9_percent_semispan_10_nodes_results_2.1_deg_AoA.csv"
    #file17_second_AoA = "94.9_percent_semispan/2.1_degrees_AoA/94.9_percent_semispan_17_nodes_results_2.1_deg_AoA.csv"
    #file25_second_AoA = "94.9_percent_semispan/2.1_degrees_AoA/94.9_percent_semispan_25_nodes_results_2.1_deg_AoA.csv"
    file45_second_AoA = "94.9_percent_semispan/2.1_degrees_AoA/94.9_percent_semispan_45_nodes_results_2.1_deg_AoA.csv"

    #File locations at an angle of attack of 4.3 degrees
    #file10_third_AoA = "94.9_percent_semispan/4.3_degrees_AoA/94.9_percent_semispan_10_nodes_results_4.3_deg_AoA.csv"
    #file17_third_AoA = "94.9_percent_semispan/4.3_degrees_AoA/94.9_percent_semispan_17_nodes_results_4.3_deg_AoA.csv"
    #file25_third_AoA = "94.9_percent_semispan/4.3_degrees_AoA/94.9_percent_semispan_25_nodes_results_4.3_deg_AoA.csv"
    file45_third_AoA = "94.9_percent_semispan/4.3_degrees_AoA/94.9_percent_semispan_45_nodes_results_4.3_deg_AoA.csv"

    #File locations at an angle of attack of 6.4 degrees
    #file10_fourth_AoA = "94.9_percent_semispan/6.4_degrees_AoA/94.9_percent_semispan_10_nodes_results_6.4_deg_AoA.csv"
    #file17_fourth_AoA = "94.9_percent_semispan/6.4_degrees_AoA/94.9_percent_semispan_17_nodes_results_6.4_deg_AoA.csv"
    #file25_fourth_AoA = "94.9_percent_semispan/6.4_degrees_AoA/94.9_percent_semispan_25_nodes_results_6.4_deg_AoA.csv"
    file45_fourth_AoA = "94.9_percent_semispan/6.4_degrees_AoA/94.9_percent_semispan_45_nodes_results_6.4_deg_AoA.csv"

    #File locations at an angle of attack of 8.6 degrees
    #file10_fifth_AoA = "94.9_percent_semispan/8.6_degrees_AoA/94.9_percent_semispan_10_nodes_results_8.6_deg_AoA.csv"
    #file17_fifth_AoA = "94.9_percent_semispan/8.6_degrees_AoA/94.9_percent_semispan_17_nodes_results_8.6_deg_AoA.csv"
    #file25_fifth_AoA = "94.9_percent_semispan/8.6_degrees_AoA/94.9_percent_semispan_25_nodes_results_8.6_deg_AoA.csv"
    file45_fifth_AoA = "94.9_percent_semispan/8.6_degrees_AoA/94.9_percent_semispan_45_nodes_results_8.6_deg_AoA.csv"

    #File locations at an angle of attack of 10.7 degrees
    #file10_sixth_AoA = "94.9_percent_semispan/10.7_degrees_AoA/94.9_percent_semispan_10_nodes_results_10.7_deg_AoA.csv"
    #file17_sixth_AoA = "94.9_percent_semispan/10.7_degrees_AoA/94.9_percent_semispan_17_nodes_results_10.7_deg_AoA.csv"
    #file25_sixth_AoA = "94.9_percent_semispan/10.7_degrees_AoA/94.9_percent_semispan_25_nodes_results_10.7_deg_AoA.csv"
    file45_sixth_AoA = "94.9_percent_semispan/10.7_degrees_AoA/94.9_percent_semispan_45_nodes_results_10.7_deg_AoA.csv"

    #Location of the experimental results
    file_experimental = "94.9_percent_semispan/94.9_percent_semispan_Experimental_Results_Wing_A.csv"
    x_LE = 46.5

    labels =["45 Nodes"]
    AoA = ["0", "2.1", "4.3", "6.4", "8.6", "10.7"]
    #labels = ["10 Nodes", "17 Nodes", "25 Nodes", "45 Nodes"]
    title = (r"Pressure Coefficient at $\frac{2y}{b} = 0.949$")


else:
    print("\n****************\nInvalid Entry. Please run script again.\n****************\n")
    quit()



#Combine file locations into an iterable list
if LE_Location == '89.8':
    file = [[file17_first_AoA, file25_first_AoA, file45_first_AoA],
        [file17_second_AoA, file25_second_AoA, file45_second_AoA],
        [file17_third_AoA, file25_third_AoA, file45_third_AoA]]

elif LE_Location == '94.9':
    file = [[file45_first_AoA],
        [file45_second_AoA],
        [file45_third_AoA],
        [file45_fourth_AoA],
        [file45_fifth_AoA],
        [file45_sixth_AoA]]

else:
    file = [[file10_first_AoA, file17_first_AoA, file25_first_AoA, file45_first_AoA],
        [file10_second_AoA, file17_second_AoA, file25_second_AoA, file45_second_AoA],
        [file10_third_AoA, file17_third_AoA, file25_third_AoA, file45_third_AoA],
        [file10_fourth_AoA, file17_fourth_AoA, file25_fourth_AoA, file45_fourth_AoA],
        [file10_fifth_AoA, file17_fifth_AoA, file25_fifth_AoA, file45_fifth_AoA],
        [file10_sixth_AoA, file17_sixth_AoA, file25_sixth_AoA, file45_sixth_AoA]]

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

        Pressure_Plot(file[i][j], labels[j], LE_Location)

    #Pull in experimental results for comparison at each angle of attack differentiating upper and lower surfaces
    upper_surface_count = Experimental[:,0].size//2 + Experimental[:,0].size%2

    plt.plot(Experimental[:upper_surface_count+1,0], Experimental[:upper_surface_count+1,i+1], ".",color="k", label="Exerimental Upper Surface")
    plt.plot(Experimental[upper_surface_count:,0], Experimental[upper_surface_count:,i+1], "*",color="k", label="Exerimental Lower Surface")



    #Plot the figure containing all curves
    plt.title(complete_title)
    plt.xlabel('x/c')
    plt.ylabel(r"$C_p$")
    plt.gca().invert_yaxis()
    if LE_Location == '94.9':
        plt.legend()
    else:
        plt.legend(ncol = 2)

    #Save the figure in its appropriate location
    # filename = LE_Location + "_percent_semispan/plots_" + LE_Location + "_percent_semispan/" + AoA[i] + "degrees_AoA_plot.pdf"
    # list_of_files.append(filename)
    # plt.savefig(filename)

    plt.show()


# Combine all pdf output files for each case in respective folder locations
# merger = PdfFileMerger()
# Combined_name = LE_Location + '_percent_semispan/plots_' + LE_Location + '_percent_semispan/Combined_plots_' + LE_Location + '_percent_semispan.pdf'
# Combined_name_copy = 'Plot_Summary/' + 'Combined_plots_' + LE_Location + '_percent_semispan.pdf'

# for k in range(len(list_of_files)):
#     merger.append(list_of_files[k])

# merger.write(Combined_name)
# merger.write(Combined_name_copy)
# merger.close()