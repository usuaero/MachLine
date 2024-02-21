from matplotlib import pyplot as plt
import numpy as np
from EURP_CURVEFITMETHOD import Angles
import math as m

# setup inputs
errorsAsPercentage = True
logTypeError = False
portNumStart = 3
portNumStart = 3
portNumEnd = 10
# aerodynamics
alpha = 20
beta = 0
static_pressure = 14.7
v_inf = 10
rho = 0.0023769

portxnum = 5
portynum = 5
portxmax = np.pi/8
portxmin = -np.pi/8
portymin, portymax = portxmin, portxmax

# setup lists

data = [[],[],[],[],[],[]]


for i in range(portNumStart,portNumEnd):
    portxnum = i
    portynum = i
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
    if errorsAsPercentage and logTypeError:
        data[0].append(i)
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
        data[0].append(i)
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
        data[0].append(i)
        data[1].append(errorList[0])
        data[2].append(errorList[1])
        data[3].append(errorList[2])
        data[4].append(errorList[3])
        data[5].append(abs(vPitot-v_inf))

xList = [i**2 for i in data[0]]
plt.plot(xList,data[1][:],label="Error Alpha")
plt.plot(xList,data[2][:],label="Error Beta")
plt.plot(xList,data[3][:],label="Error Pressure")
plt.plot(xList,data[4][:],label="Error Velocity")
plt.plot(xList,data[5][:],label="Error Pitot Velocity")
plt.xlabel("total number of ports in grid")
if errorsAsPercentage:
    plt.ylabel("% Error")
else:
    plt.ylabel("Absolute Error")
plt.title("alpha = {0} deg, beta = {1} deg".format(alpha,beta))
plt.legend()
plt.show()
    
        

