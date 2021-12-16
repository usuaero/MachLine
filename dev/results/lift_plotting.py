#This file is to take .csv or .txt files of lift distributions and plot
#them allowing comparison to experimental results

import numpy as np
import matplotlib.pyplot as plt



file = "/home/ammonhouser/gitrepository/TriPan/dev/results/Wing_A_10_nodes.csv"

Data = np.genfromtxt(file, delimiter=",", skip_header=1, dtype=float)

C_p = Data[:0]
#maximum = max(C_p)
maximum = 1
#minimum = min(C_p)
minimum = -13.5
print(C_p)
x = np.linspace(minimum-1, maximum+1, 500)

plt.figure()
plt.plot(x,C_p, label="C_p")
plt.show()