from matplotlib import pyplot as plt
import numpy as np
from EURP_CURVEFITMETHOD import Angles
import math as m

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
    errorAsPercentage = False
    aStart = -25
    aEnd = 25
    bStart = -25
    bEnd = 25
    portxnum = 5
    portynum = portxnum
    portxmax = np.radians(45)
    portxmin = -portxmax
    portymin, portymax = portxmin, portxmax

    alphas = [x for x in range(aStart,aEnd) if x != 0]
    betas = [x for x in range(bStart,bEnd) if x != 0]


    [testx, testy] = np.meshgrid(alphas,betas)
    errors = np.zeros([len(alphas),len(betas)])
    errorsP = np.zeros([len(alphas),len(betas)])

    output_file = study_dir
    output_file += "/CFM_Error_data.txt"
    with open(output_file, "w") as file:
        for i in range(len(alphas)):
            for j in range(len(betas)):
                hidden = {"rho" : 0.0023769,
                    "pinf" : 14.7,
                    "vInfMag" : 10,
                    "a":testx[0,i],
                    "b":testy[j,0]}
                input_dict = {"patern":"noCorners",
                        "numberX":portxnum,
                        "numberY":portynum,
                        "xRange":[portxmin, portxmax],
                        "yRange":[portymin, portymax],
                        "hidden": hidden
                        }
                program = Angles(input_dict)
                a, b, p, v, uncert, errorList = program.run(verbose=False,plot=False,returnErrors=True)
                
                if errorAsPercentage:
                    if a != 0:
                        errorList[0] = errorList[0]/a
                    else:
                        errorList[0] = 0
                    if b != 0:
                        errorList[1] = errorList[1]/b
                    else:
                        errorList[1] = 0
                errors[i,j] = (errorList[0]**2+errorList[1]**2)**0.5
                errorsP[i,j] = errorList[2]
                file.write(f"{errors[i,j]},")
            file.write("\n")
    file.close()



    if errorAsPercentage:
        levels = [x for x in range(1,15)] + [x for x in range(15,30,5)] + [x for x in range(30,100,10)]
    else:
        # levels = [x/5 for x in range(5)] + [x for x in range(2,10,2)]+ [x for x in range(10,100,10)]
        levels = [x/10 for x in range(10)] + [x/5 for x in range(5,10)] + [x for x in range(2,4)]
    fig, ax = plt.subplots()
    CS = ax.contour(testx,testy,errors,levels)
    if errorAsPercentage:
        ax.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=7)
    else:
        ax.clabel(CS, inline=True, fmt=fmt2, fontsize=7)

    ax.set_ylabel("Alpha (°)")
    ax.set_xlabel("Beta (°)")
    plt.show()



