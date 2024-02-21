import scipy, scipy.optimize
from math import radians
from math import factorial
import numpy as np
import json
import matplotlib as mpl
from matplotlib import pyplot as plt
import random


class Angles():
    def __init__(self,input_dict):
        import scipy, scipy.optimize
        from math import radians
        from math import factorial
        import numpy as np
        import json
        from matplotlib import pyplot as plt
        import random
        self.input_dict = input_dict
        self.getInput()

    def run(self,plot=True,verbose=True,returnErrors=False,includePitot=False):
        '''
        runs the curve fit solver to determine alpha, beta, and stagnation pressure 
        Parameters:
        plot = bool: (optional) plot ports with curve fit of pressure. Default = True
        verbose = bool: (optional) print results to screen. Default = True
        returnErrors = bool: (optional) return errors list in degrees in the form [alpha,beta,p_0]. Default = False
        Outputs:
        Alpha = np.float64: angle of attack in degrees
        Beta = np.float64: angle of sideslip in degrees (flank angle)
        P_0 = np.float64: stagnation pressure
        V = np.float64: Estimated Velocity assuming known static pressure and air density
        Uncert = list: uncertainties of each output with 1 sigma confidence
        '''
        # create lists, arrays, and counters
        ci = cj = c = 0
        self.points , c = self.generatePorts(self.patern,plot)


        # get pressure for each of the points on the sphere (sensed from sensor)
        for i in range(c):
            point = np.array([self.points[0,i],self.points[1,i],self.points[2,i]])
            pressure = self.getPressure(point)
            self.points[5,i] = pressure

        if plot:
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(projection="3d")
            # plot points
            ax2.scatter(xs=np.degrees(self.points[3,:]),ys=np.degrees(self.points[4,:]),zs=self.points[5,:],c="r",label="actual")
            ax2.set_xlabel("A port location (°)")
            ax2.set_ylabel("B port location (°)")
            ax2.set_zlabel("Pressure (psi)")
            plt.show()

        # fit curve
        initial = [1,0,0,14]
        self.params, pcov = scipy.optimize.curve_fit(self.func,[self.points[3,:],self.points[4,:]],self.points[5,:],initial,sigma=self.points[5,:])
        self.perr = np.sqrt(np.diag(pcov))
        

        self.alphaPredict = self.params[1]
        self.betaPredict = self.params[2]
        self.p_0Predict = self.params[0] + self.params[3]
        self._p_0 = 0.5 * self._rho * self._velocity**2 + self._pinf

        # this calculation assumes a known density and static pressure
        self.vPredict = (2 * (self.p_0Predict - self._pinf) / self._rho)**0.5

        # traditional pitot probe velocity
        centerPortPress = self.getPressure([1,0,0])
        vPitot = (2 * (centerPortPress - self._pinf) / self._rho)**0.5

        # find uncertainties
        self.uncert = [np.degrees(self.perr[1]),np.degrees(self.perr[2]),self.perr[0]+self.perr[3],(self.perr[0]+self.perr[3])/(2 * self.p_0Predict**0.5)]
        # find errors
        errorList = [abs(self._alpha-np.degrees(self.alphaPredict)),abs(self._beta-np.degrees(self.betaPredict)),abs(self._p_0-self.p_0Predict),abs(self._velocity-self.vPredict)]
        # print outputs
        if verbose:
            print("Parameters:\n",self.params,"\nError\n",self.perr)
            print('alpha, beta, p_0, V')
            print("Predicted")
            print(np.degrees(self.alphaPredict),np.degrees(self.betaPredict),self.p_0Predict,self.vPredict)
            print("Traditional Pitot Probe Velocity")
            print(vPitot)
            print("Uncertainty")
            print(self.uncert)
            print('input')
            print(self._alpha,self._beta,self._p_0,self._velocity)
            print('errors')
            print(errorList)
        # plot figures
        if plot:
            self.plotCurves(self.points)

        # return values
        if not returnErrors:
            return np.degrees(self.alphaPredict) , np.degrees(self.betaPredict), self.p_0Predict, self.vPredict, self.uncert
        elif includePitot:
            return self.alphaPredict, self.betaPredict, self.p_0Predict,self.vPredict, self.uncert, errorList, vPitot
        else:
            return self.alphaPredict, self.betaPredict, self.p_0Predict,self.vPredict, self.uncert, errorList

    def generatePorts(self,shape,plot):
        c = 0
        if shape == "noCorners":
            #Generate ports : number of ports
            u = v= np.linspace(np.degrees(self.rangeX[0]), np.degrees(self.rangeX[1] ) , self.numX)
            ports = []
            # generate ports : get positions of all ports 
            
            for i in u:
                for j in v:
                    if (j != 0 and i != 0):
                        if (j == i or j == -i): 
                            continue
                    x , y ,z = self.alphaBetaToCart(np.radians(i), np.radians(j))
                    port = [x,y,z,np.radians(i),np.radians(j),0]
                    ports.append(port)
                    c += 1
            # add transposed dictionary to ports
            points = np.transpose(np.array(ports))


        elif shape == "triad":
            ci = 0
            ports = []
            for i in range(self.numX):
                '''generates an orthogonal triad of port positions on a unit sphere based on 3 euler angles pitch, roll, and yaw in degrees'''
                # inputs
                
                pitchAngle = np.radians(pitchAngle)  #fix this to take as input
                rollAngle = np.radians(rollAngle)
                yawAngle = np.radians(yawAngle)
                # generate transformation matricies
                yaw = np.array([[np.cos(yawAngle),-np.sin(yawAngle),0],
                                [np.sin(yawAngle),np.cos(yawAngle),0],
                                [0,0,1]])
                pitch = np.array([[np.cos(pitchAngle),0,np.sin(pitchAngle)],
                                [0,1,0],
                                [-np.sin(pitchAngle),0,np.cos(pitchAngle)]])
                roll = np.array([[1,0,0],
                                [0,np.cos(rollAngle),-np.sin(rollAngle)],
                                [0,np.sin(rollAngle),np.cos(rollAngle)]])
                # calculate full transformation matrix 
                t = yaw @ pitch @ roll
                # get ports
                p = t @ np.eye(3) 
                # save matrix to allow for reorganizing
                self.tfs.append(p)
                # add ports to ports dictionary
                for i in range(3):
                    a , b = self.cartToAlphaBeta(p[0,i],p[1,i],p[2,i])
                    self.port = [p[0,i],p[1,i],p[2,i],a,b]
                    ports.append(port)
                    ci += 1
                
            
        else:
            if shape != "grid":
                print("Error in Port Shape Generation, using default grid shape")
            #Generate ports : number of ports
            u = np.linspace(self.rangeX[0], self.rangeX[1]  , self.numX)
            v = np.linspace(self.rangeY[0], self.rangeY[1] , self.numY)
            comb = self.numX * self.numY
            points = np.zeros([6,comb])

            # generate ports : get positions of all ports 
            for i in u:
                for j in v:
                    x , y ,z = self.alphaBetaToCart(i, j)
                    points[0,c] = x
                    points[1,c] = y
                    points[2,c] = z
                    points[3,c] = i
                    points[4,c] = j
                    c += 1
        
        if plot:
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(projection="3d")
            # plot points
            ax1.scatter(xs=points[0,:],ys=points[1,:],zs=points[3,:],c="r",label="Ports")
            plt.show()

        return points, c
        


    def getInput(self):
        # set tolerance for rounding
        self.tol = 10
        # get range for grid patern
        self.patern = self.input_dict["patern"]
        if self.patern == "grid":
            self.numX = self.input_dict['numberX']
            self.numY = self.input_dict['numberY']
            self.rangeX = self.input_dict["xRange"]
            self.rangeY = self.input_dict["yRange"]
        elif self.patern == "noCorners":
            self.numX = self.input_dict['numberX']
            self.rangeX = self.input_dict["xRange"]
            self.rangeY = self.rangeX
        
        
        # input hidden values
        self._velocity = self.input_dict["hidden"]["vInfMag"]
        self._alpha = self.input_dict["hidden"]["a"]
        self._beta = self.input_dict["hidden"]["b"]
        self._rho = self.input_dict["hidden"]["rho"]
        self._pinf = self.input_dict["hidden"]["pinf"]
        # Calculate Hidden Values
        vx , vy, vz = self.alphaBetaToCart(np.radians(self._alpha), np.radians(self._beta))
        self._vinf = [vx, vy, vz]
        self._p_0 = 0.5 * self._rho * self._velocity**2 + self._pinf


    def alphaBetaToCart(self,a,b):
        '''
        Outputs Cartesian values assuming a radius of 1 unit
        Parameters:
        a = float: alpha location in radians.
        b = float: beta location in radians.
        Outputs: 
        x = np.float64: cartesian x location.
        y = np.float64: cartesian y location.
        z = np.float64: cartesian z location.
        '''
        x = np.cos(a)*np.cos(b)/(1-np.sin(a)**2*np.sin(b)**2)**0.5
        y = np.cos(a)*np.sin(b)/(1-np.sin(a)**2*np.sin(b)**2)**0.5
        z = np.sin(a)*np.cos(b)/(1-np.sin(a)**2*np.sin(b)**2)**0.5
        return x,y,z


    def cartToAlphaBeta(self,x,y,z):
        '''
        Outputs Alpha and Beta values assuming a radius of 1 unit
        Parameters:
        x = float: cartesian x location. 
        y = float: cartesian y location. 
        z = float: cartesian z location. 
        Outputs: 
        a = np.float64: alpha location in radians. 
        b = np.float64: beta location in radians. 
        '''
        a = np.arctan2(z,x)
        b = np.arctan2(y,x)
        return a, b

    def func(self,coords,a,b,c,d):
        x = coords[0]
        y = coords[1]
        ans = a * np.cos(((x-b)**2+(y-c)**2)**0.5) + d
        return ans

    # def func(self,coords,a,b,c,d):
    #     x = coords[0]
    #     y = coords[1]
    #     ans = a*((x-b)**2+(y-c)**2) + d
    #     return ans

    def getPressure(self, point):
        '''
        Outputs the pressure of a point using the hidden values inported on init
        Parameters: 
        point = list:  Cartisian port location.
        Outputs = float: pressure 
        '''

        mag = np.dot(self._vinf,point)
        t = np.arccos(round(mag,self.tol))
        cp = 1-9/4*(np.sin(t))**2
        pressure = cp * 0.5 * self._rho * self._velocity**2 + self._pinf


        return pressure
        

       
    def plotCurves(self,points):
        # plotting section.  Will be moved to seperate method soon
        # get test data to plot
        ci = cj = c = 0
        maxx = 50
        maxy = 50
        u = np.linspace(self.rangeX[0], self.rangeX[1] , maxx)
        v = np.linspace(self.rangeY[0], self.rangeY[1], maxy)
        fit = np.zeros([3,maxx*maxy])

        # get points on sphere
        for i in u:
            for j in v:
                fit[0,c] = i
                fit[1,c] = j
                fit[2,c] = self.func([i,j], self.params[0], self.params[1], self.params[2], self.params[3])
                c += 1

        xp, yp = fit[0,:], fit[1,:]
        zp = fit[2,:]
        # create figure
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")

        # plot points
        actPlot = ax.scatter(xs=np.degrees(points[3,:]),ys=np.degrees(points[4,:]),zs=points[5,:],c="r",label="actual")

        # plot surface from curve fit
        fitPlot = ax.plot_trisurf(np.degrees(xp),np.degrees(yp),zp,antialiased=True,label="curve fit")
        # self.ax3.set_box_aspect((np.ptp(points[0,:]), np.ptp(points[1,:]), np.ptp(points[2,:])))
        ax.set_xlabel("A port location (°)")
        ax.set_ylabel("B port location (°)")
        ax.set_zlabel("Pressure (psi)")
        # fixing a plotting error pyplot
        fitPlot._edgecolors2d = fitPlot._edgecolor3d
        fitPlot._facecolors2d = fitPlot._facecolor3d
        plt.legend()
        plt.show()

if __name__ == "__main__":
    mpl.rcParams['figure.dpi'] = 200
    singleRun = True
    aoa2d = False
    ab3d = False
    contour = False
    # run once
    portxnum = 3
    portxmax = radians(45)
    portxmin = -portxmax
    portynum , portymin , portymax = portxnum, portxmin , portxmax
    if singleRun:
        hidden = {"rho" :0.0023769,
            "pinf" : 14.7,
            "vInfMag" : 10,
            "a":10,
            "b":5}
        input_dict = {"patern":"noCorners",
                "numberX":portxnum,
                "numberY":portynum,
                "xRange":[portxmin, portxmax],
                "yRange":[portymin, portymax],
                "hidden": hidden
                }
        program = Angles(input_dict)
        a, b, p, v, uncert, errorList = program.run(verbose=True,plot=True,returnErrors=True)

    # 2d plot error vs aoa
    if aoa2d:
        alphas = []
        errorsA = []
        errorsB = []
        errorsP = []
        for i in range(20):
            hidden = {"rho" : 0.0023769,
                "pinf" : 14.7,
                "vInfMag" : 10,
                "a":i,
                "b":0}
            input_dict = {"patern":"grid",
                    "numberX":5,
                    "numberY":5,
                    "xRange":[-np.pi/8, np.pi/8],
                    "yRange":[-np.pi/8, np.pi/8],
                    "hidden": hidden
                    }
            program = Angles(input_dict)
            a, b, p,v, uncert, errorList = program.run(verbose=False,plot=False,returnErrors=True)
            alphas.append(i)
            errorsA.append(errorList[0])
            errorsB.append(errorList[1])
            errorsP.append(errorList[2])
            
        
        plt.plot(alphas,errorsA,label="Alpha Error")
        plt.plot(alphas,errorsB,label="Beta Error")
        plt.plot(alphas,errorsP,label="Pressure Error")
        plt.legend()
        plt.show()
    # 3d plot with alpha and beta 
    if ab3d:
        errorAsPercentage = True
        aStart = -15
        aEnd = 15
        bStart = -15
        bEnd = 15
        portxnum = 3
        portynum = 3
        portxmax = np.pi/4
        portxmin = -np.pi/4
        portymin, portymax = portxmin, portxmax

        alphas = np.linspace(aStart,aEnd,aEnd-aStart+1)
        betas = np.linspace(bStart,bEnd,bEnd-bStart+1)
        alist = []
        blist = []
        errorsA = []
        errorsB = []
        errorsP = []
        for i in range(len(alphas)):
            for j in range(len(betas)):
                hidden = {"rho" : 0.000001375520833,
                    "pinf" : 14.7,
                    "vInfMag" : 10,
                    "a":alphas[i],
                    "b":betas[j]}
                input_dict = {"patern":"grid",
                        "numberX":portxnum,
                        "numberY":portynum,
                        "xRange":[portxmin, portxmax],
                        "yRange":[portymin, portymax],
                        "hidden": hidden
                        }
                program = Angles(input_dict)
                a, b, p, v, uncert, errorList = program.run(verbose=False,plot=False,returnErrors=True)
                alist.append(alphas[i])
                blist.append(betas[j])
                if errorAsPercentage:
                    if a != 0:
                        errorsA.append(abs(errorList[0]/a))
                    else: 
                        errorsA.append(0)
                    if b != 0:
                        errorsB.append(abs(errorList[1]/b))
                    else:
                        errorsB.append(0)
                    if p != 0:
                        errorsP.append(abs(errorList[2]/p))
                    else:
                        errorsP.append(0)
                else:
                    errorsA.append(errorList[0])
                    errorsB.append(errorList[1])
                    errorsP.append(1000*errorList[2])
            
        
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        ax.scatter(alist,blist,errorsA,label="Alpha Error (°)")
        ax.scatter(alist,blist,errorsB,label="Beta Error (°)")
        ax.scatter(alist,blist,errorsP,label="Pressure Error (0.001psi)")
        ax.set_xlabel("Alpha")
        ax.set_ylabel("Beta")
        ax.set_zlabel("Error")
        plt.legend()
        plt.show()
    
    if contour:
        errorAsPercentage = False
        aStart = -15
        aEnd = 15
        bStart = -15
        bEnd = 15
        portxnum = 5
        portynum = 5
        portxmax = np.pi/8
        portxmin = -np.pi/8
        portymin, portymax = portxmin, portxmax

        alphas = np.linspace(aStart,aEnd,aEnd-aStart+1)
        betas = np.linspace(bStart,bEnd,bEnd-bStart+1)
        
        [testx, testy] = np.meshgrid(alphas,betas)
        errors = np.zeros([len(alphas),len(betas)])
        errorsP = np.zeros([len(alphas),len(betas)])
        
        for i in range(len(alphas)):
            for j in range(len(betas)):
                hidden = {"rho" : 0.0023769,
                    "pinf" : 14.7,
                    "vInfMag" : 10,
                    "a":testx[0,i],
                    "b":testy[j,0]}
                input_dict = {"patern":"grid",
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
        
        levels = np.linspace(np.amin(errors), np.amax(errors),10)
        plt.contour(testx,testy,errors,levels=20)
        plt.colorbar()
        plt.show()







    





            