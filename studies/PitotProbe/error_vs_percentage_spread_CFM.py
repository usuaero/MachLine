from matplotlib import pyplot as plt
import numpy as np
from EURP_CURVEFITMETHOD import Angles
import math as m

# setup inputs
errorsAsPercentage = True
percentAlpha = False
logTypeError = False
portNumStart = 3
portNumEnd = 10
portSpreadStart = 15
portSpreadEnd = 60
degPerPort = 15


# aerodynamics
static_pressure = 14.7
v_inf = 10
rho = 0.0023769
beta=0

# setup lists

data = [[],[],[],[],[],[]]
bigData = []
names = []

for i in range(portSpreadStart,portSpreadEnd+1,5):
    portxnum = int(2 * i / degPerPort) + 1
    portynum = portxnum
    for j in range(1,45):
        alpha = j
        
        portxmax = np.radians(i)
        portxmin = -portxmax
        portymin, portymax = portxmin, portxmax
        dynamic_pressure = 1/2 * rho * v_inf **2
        hidden = {"rho" : rho,
            "pinf" : static_pressure,
            "vInfMag" : 10,
            "a":alpha,
            "b":beta}
        input_dict = {"patern":"grid",
                "numberX":portxnum,
                "numberY":portynum,
                "xRange":[portxmin, portxmax],
                "yRange":[portymin, portymax],
                "hidden": hidden
                }
        program = Angles(input_dict)
        a, b, p, v, uncert, errorList, vPitot = program.run(verbose=False,plot=False,returnErrors=True,includePitot=True)
        if percentAlpha:
            data[0].append(j/i*100)
        else:
            data[0].append(j)

        if errorsAsPercentage and logTypeError:
            alphaError = m.log(m.e**(alpha - errorList[0]))
            betaError = m.log(m.e**(beta - errorList[1]))
            pressureError = m.log(m.e**(dynamic_pressure - errorList[2]))
            velocityError = m.log(m.e**(v_inf - errorList[3]))
            pitotVelocityError = m.log(m.e**(v_inf - abs(vPitot-v_inf)))

            data[1].append(alphaError)
            data[2].append(betaError)
            data[3].append(pressureError)
            data[4].append(velocityError)
            data[5].append(pitotVelocityError)
        elif errorsAsPercentage:
            if alpha == 0:
                alphaError = 0
            else:
                alphaError = errorList[0]/alpha *100
            if beta == 0:
                betaError = 0
            else:    
                betaError = errorList[1]/beta * 100

            pressureError = errorList[2]/dynamic_pressure * 100
            velocityError = errorList[3]/v_inf * 100
            pitotVelocityError = abs(vPitot-v_inf)/v_inf * 100

            data[1].append(alphaError)
            data[2].append(betaError)
            data[3].append(pressureError)
            data[4].append(velocityError)
            data[5].append(pitotVelocityError)
        else:
            data[1].append(errorList[0])
            data[2].append(errorList[1])
            data[3].append(errorList[2])
            data[4].append(errorList[3])
            data[5].append(abs(vPitot-v_inf))
    bigData.append(data)
    data = [[],[],[],[],[],[]]
    names.append([i,portxnum])


for i in range(len(bigData)):
    xList = bigData[i][0]
    plt.scatter(xList,bigData[i][1],label="{0} deg spread, {1} count".format(names[i][0],names[i][1]))
if percentAlpha:
    plt.xlabel("Alpha (% of spread)")
else:
    plt.xlabel("Alpha (deg)")
if errorsAsPercentage:
    plt.ylabel("% Alpha Error")
else:
    plt.ylabel("Absolute Alpha Error")
plt.legend()
plt.show()

for i in range(len(bigData)):
    xList = bigData[i][0]
    plt.plot(xList,bigData[i][4],label="{0} deg spread, {1} count".format(names[i][0],names[i][1]))
if percentAlpha:
    plt.xlabel("Alpha (% of spread)")
else:
    plt.xlabel("Alpha (deg)")
if errorsAsPercentage:
    plt.ylabel("% Velocity Error")
else:
    plt.ylabel("Absolute Velocity Error")
plt.legend()
plt.show()
        
            

