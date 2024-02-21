from matplotlib import pyplot as plt
import numpy as np
from EURP_PEM import Sphere
import math as m
import csv

def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"

def fmt2(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \°" if plt.rcParams["text.usetex"] else f"{s} °"

if __name__ == "__main__":
    study_dir = "studies/PitotProbe"
    run_sim= False
    errorAsPercentage = True
    aStart = -25
    aEnd = 25
    bStart = -25
    bEnd = 25
    portxnum = 5
    portynum = portxnum
    portxmax = np.radians(45)
    portxmin = -portxmax
    portymin, portymax = portxmin, portxmax

    alphas = [x for x in range(aStart,aEnd) ]
    betas = [x for x in range(bStart,bEnd) ]


    [testx, testy] = np.meshgrid(alphas,betas)
    errors = np.zeros([len(alphas),len(betas)])
    
    output_file = study_dir
    output_file += "/PEM_Error_data.txt"
    if run_sim:
        
        with open(output_file, "w") as file:
            for i in range(len(alphas)):
                for j in range(len(betas)):
                    
                    study_directory = "studies/PitotProbe"
                    program = Sphere(study_dir,"input_EURP.json",alpha=alphas[i],beta=betas[j])
                    a, b, aError, bError = program.run()
                    # aError = i
                    # bError = j
                    errorMag = (aError**2+bError**2)**0.5
                    errors[i,j] = errorMag
                    file.write(f"{errorMag},")
                file.write("\n")
        file.close()


    with open(output_file, 'r') as csv_file:
        reader = csv.reader(csv_file)
        num_lines = 50
        print(num_lines)
        row_number = 0
        for row in reader:            
            for i in range(len(row)):
                if errorAsPercentage:
                    if ((alphas[row_number]**2+betas[i]**2)**0.5) == 0:
                        errors[row_number,i] = 0
                    else:
                        errors[row_number,i] = float(row[i])/((np.radians(alphas[row_number])**2+np.radians(betas[i])**2)**0.5)*100
                        
                else:
                    errors[row_number,i] = float(row[i]) 
            row_number +=1
    csv_file.close()
    if errorAsPercentage:
        levels = [x for x in range(1,15)]
    else:
        #levels = [x/5 for x in range(5)] + [x for x in range(2,10,2)]+ [x for x in range(10,100,10)]
        levels = [x/10 for x in range(10)] + [x/5 for x in range(5,10)] + [x for x in range(2,4)]
        errors = np.degrees(errors)
    fig, ax = plt.subplots()
    CS = ax.contour(testx,testy,errors,levels)
    if errorAsPercentage:
        ax.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=7)
    else:
        ax.clabel(CS, inline=True, fmt=fmt2, fontsize=7)

    ax.set_ylabel("Alpha (°)")
    ax.set_xlabel("Beta (°)")
    plt.show()



