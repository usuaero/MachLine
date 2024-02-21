from matplotlib import pyplot as plt
import numpy as np
from EURP_CURVEFITMETHOD import Angles
import math as m

# setup inputs
errorsAsPercentage = True
logTypeError = False
portNumStart = 3
portNumEnd = 10
portSpreadStart = 5
portSpreadEnd = 85
# aerodynamics
alpha = 10
beta = 5
static_pressure = 14.7
v_inf = 10
rho = 0.0023769

# setup lists

count = range(portSpreadStart,portSpreadEnd,2)
spread = range(portNumStart,portNumEnd)
data=[[],[],[]]
[testx, testy] = np.meshgrid(spread,count)
errorsA = np.zeros([len(count),len(spread)])
errorsB = np.zeros([len(count),len(spread)])
errorsP = np.zeros([len(count),len(spread)])
errorsV = np.zeros([len(count),len(spread)])
density = np.zeros([len(count),len(spread)])

for i in range(len(count)):
    for j in range(len(spread)):
        portxnum = testx[i,j]
        portynum = portxnum
        portxmax = np.radians(testy[i,j])
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
        data[0].append(count[i])
        data[1].append(spread[j])
        data[2].append(errorList[0])
        arcArea = np.radians(2* portxmax) * np.radians(2 *  portymax)
        portDensity = portxnum**2/arcArea
        density[i,j] = portDensity
        if errorsAsPercentage and logTypeError:
            alphaError = m.log(m.e**(alpha - errorList[0]))
            betaError = m.log(m.e**(beta - errorList[1]))
            pressureError = m.log(m.e**(dynamic_pressure - errorList[2]))
            velocityError = m.log(m.e**(v_inf - errorList[3]))
            pitotVelocityError = m.log(m.e**(v_inf - abs(vPitot-v_inf)))

            errorsA[i,j] = alphaError
            errorsB[i,j] = betaError
            errorsP[i,j] = pressureError
            errorsV[i,j] = velocityError
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

            errorsA[i,j] = alphaError
            errorsB[i,j] = betaError
            errorsP[i,j] = pressureError
            errorsV[i,j] = velocityError
            
        else:
            alphaError = errorList[0]
            betaError = errorList[1]
            pressureError = errorList[2]
            velocityError = errorList[3]

            errorsA[i,j] = alphaError
            errorsB[i,j] = betaError
            errorsP[i,j] = pressureError
            errorsV[i,j] = velocityError
testy = testy*2
plt.title("Port Density")
# alpha error
CS = plt.contour(testx,testy,density,levels=50)
plt.clabel(CS, inline=True, fontsize=7)
plt.xlabel("Horizontal Number of Ports")
plt.ylabel("Spread Angle(deg)")
plt.show()

plt.title("Alpha error for alpha = {0} deg, beta = {1} deg".format(alpha,beta))
# alpha error
CS = plt.contour(testx,testy,errorsA,levels=40)
plt.clabel(CS, inline=True, fontsize=7)
plt.xlabel("Horizontal Number of Ports")
plt.ylabel("Spread Angle(deg)")
plt.show()

# beta error
plt.title("Beta error for alpha = {0} deg, beta = {1} deg".format(alpha,beta))
CS = plt.contour(testx,testy,errorsB,levels=50)
plt.clabel(CS, inline=True, fontsize=7)
plt.xlabel("Horizontal Number of Ports")
plt.ylabel("Spread Angle(deg)")
plt.show()

# velocity  error
plt.title("Velocity error for alpha = {0} deg, beta = {1} deg".format(alpha,beta))
CS = plt.contour(testx,testy,errorsV,levels=30)
plt.clabel(CS, inline=True, fontsize=10)
plt.xlabel("Horizontal Number of Ports")
plt.ylabel("Spread Angle(deg)")
plt.show()



# plt.scatter(data[0],data[1],c=data[2])
# plt.colorbar()
# plt.show()
    
        

