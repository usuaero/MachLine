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
    s = f"{x:.2f}"
    if s.endswith("0"):
        s = f"{x:.1f}"
    return rf"{s} \°" if plt.rcParams["text.usetex"] else f"{s} °"

if __name__ == "__main__":
    study_dir = "studies/PitotProbe"
    run_sim= True
    errorAsPercentage = True
    

    
    alphas = [x for x in range(-25,25) ]
    betas = [x for x in range(-25,25) ]
    [testx, testy] = np.meshgrid(alphas,betas)
    sigmaA = np.zeros([len(alphas),len(betas)])
    sigmaB = np.zeros([len(alphas),len(betas)])
    output_file = study_dir
    output_file += "/PEM_Error_data_monteCarlo.txt"
    if run_sim:
        with open(output_file, "w") as file:
            for i in range(len(alphas)):
                for j in range(len(betas)):
                    errorsA = []
                    errorsB = []
                    for k in range(50):
                        program = Sphere(study_dir,"input_EURP.json",alpha=alphas[i],beta=betas[j])
                        a, b, aError, bError = program.run()
                        # aError = i
                        # bError = j
                        errorsA.append(np.degrees(aError))
                        errorsB.append(np.degrees(bError))
                        # file.write(f"{a},{b},{aError},{bError}\n")
                    sigmaA[i,j] = np.std(errorsA)
                    sigmaB[i,j] = np.std(errorsB)
                    file.write(f"{sigmaA[i,j]},{sigmaB[i,j]}\n")
        file.close()

    levels = [x/20 for x in range(10)]
    fig, ax = plt.subplots()
    CS = ax.contour(testx,testy,sigmaA,levels)
    ax.clabel(CS, inline=True, fmt=fmt2, fontsize=7)
    # ax.scatter(errorsB,errorsA)
    ax.set_ylabel("Alpha (°)")
    ax.set_xlabel("Beta (°)")
    plt.show()


    fig2, ax2 = plt.subplots()
    CS = ax2.contour(testx,testy,sigmaB,levels)
    ax2.clabel(CS, inline=True, fmt=fmt2, fontsize=7)
    # ax.scatter(errorsB,errorsA)
    ax2.set_ylabel("Alpha (°)")
    ax2.set_xlabel("Beta (°)")
    plt.show()



