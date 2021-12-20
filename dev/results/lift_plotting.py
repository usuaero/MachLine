#This file is to take .csv or .txt files of lift distributions and plot
#them allowing comparison to experimental results

import numpy as np
import matplotlib.pyplot as plt

def Pressure_Plot(file, LABEL):
    Data = np.genfromtxt(file, delimiter=",", skip_header=1, dtype=float)

    C_p = -1* Data[:,0]
    x = Data[:,1]
    
    plt.plot((x-24.5)/20,C_p, label=LABEL)
    plt.xlabel('x/c')
    plt.ylabel(r"$-C_p$")
    plt.legend()
    plt.title("Pressure Coefficient at Mid Semispan")

    return

#Main
#Input Files
file10 = "/home/ammonhouser/gitrepository/TriPan/dev/results/Wing_A_10_nodes.csv"
file15 = "/home/ammonhouser/gitrepository/TriPan/dev/results/Wing_A_15_nodes.csv"
file20 = "/home/ammonhouser/gitrepository/TriPan/dev/results/Wing_A_20_nodes.csv"
file25 = "/home/ammonhouser/gitrepository/TriPan/dev/results/Wing_A_25_nodes.csv"
file45 = "/home/ammonhouser/gitrepository/TriPan/dev/results/Wing_A_45_nodes.csv"

file = [file10, file15, file20, file25, file45]

#Labels for curves on the plot
labels = ["10 Nodes", "15 Nodes", "20 Nodes", "25 Nodes", "45 Nodes"]

#Call Pressure_Plot for all Node files
for i in range(len(file)):

    Pressure_Plot(file[i], labels[i])

plt.show()


# print("before quit")
# breakpoint()
# quit()
# print("after quit")