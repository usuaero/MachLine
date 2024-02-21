# Hurwitz EURP alpha and beta sim
#imports
from math import radians
import numpy as np
import json
import matplotlib as mpl
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
import random
import scipy.optimize
import csv
import get_pressure_Dists as ML

class Sphere():
    def __init__(self,study_directory,input_filename,alpha=0,beta=0):
        self.debug = False
        self.study_directory = study_directory
        self.filename = study_directory + "/" + input_filename
        self.alpha_in = alpha
        self.beta_in = beta
        # this is the number of decimal places rounded off to midigate numerical trig error
        self.tol = 10
        # initialize hidden variables used for verificaiton 
        self._thetas = []
        self._pressures = []
        self._cps = []
        self._vinf = []

    
            

    
    
    def getAB(self):
        '''This is the main function used to determine alpha and beta.
        This function has no inputs and is run using the run function.  The run function handles simulation and graphing.
        '''
        # first get pressure from ports and save to ports dictionary and pressures list
        # calculate the pressures at each port from machline
        if self.machline:
            self._getMachlinePressure()
            
        else:
        # get pressures from ports
            for i in self.portsKeys:
                p = self._getP(i)
                self.ports[i]["P"] = p
                self.pressures.append(p)
        # method to get cp from pressures to go here for now actual cp is used

        # __________code to add actual cp to ports dictionary_______________
        # get stagnation pressure (by cheating)
        self.actualp_0 = 0.5 * self.rho * self._vinfmag**2 + self.pinf
        self.p_0 = self.actualp_0
        # self.p_0, self.fitAlpha, self.fitBeta = self.getP0()
        # calculate cp from given stagnation pressure (this bit is not cheating, actual code)
        for i in self.portsKeys:
            index = list(self.portsKeys).index(i)
            cp = (self.ports[i]["P"] - self.pinf) / (self.p_0 - self.pinf)
            self.ports[i]["CP"] = cp
            self.cps.append(cp)
        # __________end code to add actual cp to ports dicationary__________
        
        # find theta for each port 
        for i in self.portsKeys:
            cp = self.ports[i]["CP"]
            t = np.arcsin((4/9*(1-cp))**0.5)
            self.ports[i]["theta"] = t
            self.thetas.append(t)

        # use given triads to get v, alpha, and beta in each frame
        vpredicttriad = np.zeros(3)
        self.vpredict = []
        self.alphaPredict = []
        self.betaPredict = []
        for i in range(len(self.triads)):
            
            for j in range(3):
                name = "point" + str(i) + str(j+1)
                vpredicttriad[j] = np.cos(self.ports[name]['theta'])
            self.vpredict.append(self.tfs[i] @ vpredicttriad)
            pAlpha, pBeta = self.cartToAlphaBeta(self.vpredict[i][0], self.vpredict[i][1], self.vpredict[i][2])
            self.alphaPredict.append(pAlpha)
            self.betaPredict.append(pBeta)


    def func(self,coords,a,b,c,d):
        x = coords[0]
        y = coords[1]
        ans = a * np.cos(((x-b)**2+(y-c)**2)**0.5) + d
        return ans

    def getP0(self):
        initial = [1,0,0,14]
        alphas = []
        betas = []
        ps = []
        for i in self.portsKeys:
            pos = self.ports[i]['position']
            p = self.ports[i]['P']
            a,b = self.cartToAlphaBeta(pos[0], pos[1], pos[2])
            if a <= np.pi/4 and b <= np.pi/4:
                alphas.append(a)
                betas.append(b)
                ps.append(p)
        self.params, pcov = scipy.optimize.curve_fit(self.func,[alphas,betas],ps,initial)
        perr = np.sqrt(np.diag(pcov))
        p_0 = self.params[0] + self.params[3]
        alpha = self.params[1]
        beta = self.params[2]
        return p_0 , alpha , beta
        



    #  ___________________________  EVERYTHING BELOW THIS LINE IS "HIDDEN" AND IS PART OF THE SIMULATION, NOT THE ALPHA AND BETA METHOD____________________
    
    
    def run(self):
        # get inputs
        self.getInput()
        # create lists
        self.pressures = []
        self.cps = []
        self.thetas = []

        
        # find alpha and beta
        self.getAB()



        self.pressures = np.array(self.pressures)
        self.thetas = np.array(self.thetas)

        

        # plot sphere and ports figure
        if self.plotPorts and not self.plotSphere:
            self.fig2 = plt.figure(1)
            self.ax2 = plt.subplot(projection='3d')
            self._plotPorts()
            self.ax2.set_xlabel("x-axis")
            self.ax2.set_ylabel("y-axis")
            self.ax2.set_zlabel("z-axis")
            

        if self.plotSphere:
           self._plotSphere()
           self._plotPorts()
           
        if self.arrows:
            vx,vy,vz = zip(self._vinf)
            plt.quiver(0,0,0,vx,vy,vz,color = "red", label="Actual")
            for i in range(len(self.triads)):
                name = "Triad " + str(i)
                vx,vy,vz = zip(self.vpredict[i])
                plt.quiver(0,0,0,vx,vy,vz,label=name)
        if self.arrows or self.plotSphere or self.plotPorts:
            plt.legend()
            plt.show()

        if self.mapping:
            self._plotMapping()


        # plot sine wave figure with thetas for ports
        if self.plotCPvTheta:
            fig1 = plt.figure(2)
            ax1 = plt.subplot()
            ax1.scatter(x=np.degrees(self.thetas),y=self.cps,label="points")
            thetaIdeal = np.linspace(0,np.pi,100)
            pressureIdeal = 1-9/4*(np.sin(thetaIdeal))**2
            ax1.plot(np.degrees(thetaIdeal),pressureIdeal,color="orange",label="Theory")
            ax1.set_xlabel("Theta")
            ax1.set_ylabel("Cp")
            plt.show()
        


        if self.debug:
            print("pressures\n", self.pressures)
            print("actual pressures\n", self._pressures)
            print("\nPredicted P_0\n",self.p_0)
            print("Actual P_0\n",self.actualp_0)
            print("\ncp\n", self.cps)
            print("actual cp\n", self._cps)
            print("\ntheta\n", self.thetas)
            print("actual theta\n", self._thetas)
            print("\npredicted vinfs\n", self.vpredict)
            print("actual vinf\n", self._vinf)
            print("\nactual alpha\n",np.degrees(self._alpha))
            print("predicted alpha\n", np.degrees(self.alphaPredict))
            print("averaged predicted alpha\n", np.nanmean(np.degrees(self.alphaPredict)))
            # print("Curve fit Alpha\n",np.degrees(self.fitAlpha))
            print("\nactual beta\n", np.degrees(self._beta))
            print("predicted beta\n", np.degrees(self.betaPredict))
            print("averaged predicted beta\n", np.nanmean(np.degrees(self.betaPredict)))
            # print("Curve fit beta\n",np.degrees(self.fitBeta))

            print("Port Locations\n")
            for i in self.portsKeys:
                print(self.ports[i]['position'])

        errorA = self._alpha - np.nanmean(self.alphaPredict)
        errorB = self._beta - np.nanmean(self.betaPredict)
        return np.nanmean(self.alphaPredict), np.nanmean(self.betaPredict), errorA, errorB


    def getInput(self):
        # get file
        self.json_string = open(self.filename).read()
        self.input_dict = json.loads(self.json_string)

        # import parameters
        self.radius = self.input_dict["sphere"]["radius[in]"]
        self.machline = self.input_dict["settings"]["machline"]
        self.noise = self.input_dict["settings"]["noise"]

        # get port locations
        self.triads = list(self.input_dict['triads'].keys())
        self.ports = dict()
        self.tfs = []
        for i in self.triads:
            p = self.input_dict['triads'][i]['pitch']
            r = self.input_dict['triads'][i]['roll']
            y = self.input_dict['triads'][i]['yaw']
            self.genPorts(i, p, r, y)
        self.portsKeys = self.ports.keys()

        # get aerodynamic properties (hidden)
        self._vinfmag = self.input_dict["freestream"]["velocity[ft/s]"]
        alphaDeg = self.input_dict["freestream"].get("alpha[deg]",self.alpha_in)
        betaDeg = self.input_dict["freestream"].get("beta[deg]",self.beta_in)
        self.rho = self.input_dict["freestream"].get("density[slugs/ft^3]",0.0023769)
        self.pinf = self.input_dict["freestream"].get("pressure[psi]",14.7)

        # import plot options
        self.plotSphere = self.input_dict["plot"]["spherical_pressure"]
        self.plotPorts = self.input_dict["plot"]["ports"]
        self.plotCPvTheta = self.input_dict["plot"]["CPvTheta"]
        self.arrows = self.input_dict["plot"]["arrows"]
        self.mapping = self.input_dict["plot"]["pressure_mapping"]

        # Switch units to radians
        self._alpha = radians(alphaDeg)
        self._beta = radians(betaDeg)


    def alphaBetaToCart(self,a,b):
        x = np.cos(a)*np.cos(b)
        y = np.sin(b)
        z = np.sin(a)*np.cos(b)
        return x,y,z
    
    def cartToAlphaBeta(self,x,y,z):
        a = np.arctan(z/x)
        b = np.arcsin(y/(x**2+y**2+z**2)**0.5)
        return a, b

    def _getTheta(self,portName):
            '''this function determines the actual theta value using the given actual alpha and beta and returns the theta value from a given port name'''
            # get port position
            
            portDirection = self.ports[portName]["position"]
            # get velocity direction
            vx , vy, vz = self.alphaBetaToCart(self._alpha,self._beta)
            vDirection = np.array([vx,vy,vz])
            self._vinf = vDirection
            mag = np.dot(vDirection,portDirection)
            t = np.arccos(round(mag,self.tol))
            return t
        
    def _getP(self,portName):
        '''This function determines the actual C_p and pressure for a given port using the _getTheta Function and the vinfmag value from the input.  It then returns the pressure indicated at the port with added noise'''
        # get actual theta for given port
        theta = self._getTheta(portName)
        self._thetas.append(theta)

        # get actual cp for port
        cp = 1-9/4*(np.sin(theta))**2
        
        self._cps.append(cp)

        # get actual pressure for port
        pressure = cp * 0.5 * self.rho * self._vinfmag**2 + self.pinf
        self._pressures.append(pressure)

        # addition of pressure sensor noise 
        if self.noise:
             pressure += random.randint(-10,10)/10 * 0.00050763
            #pressure += 0.00050763
        return pressure

    def _calcMachlinePressure(self):
        # export port locations for machline run
        self._exportPorts()
        # run machline to get pressures
        N_sys, l_avg, C_F, C_M = ML.run_machline_for_ab(self._alpha,self._beta,self.study_directory)


    def _getMachlinePressure(self):
        self._calcMachlinePressure()
        x,y,z,CP = self._importPortPressures()
        pressures = []
        index = 0
        for i in self.portsKeys:
            theta = self._getTheta(i)
            self._thetas.append(theta)
            for j in range(len(x)):
                if (round(x[j],6) == round(self.ports[i]["position"][0],6) and round(y[j],6) == round(self.ports[i]["position"][1],6) and round(z[j],6) == round(self.ports[i]["position"][2],6)):
                    pressure = CP[j] * 0.5 * self.rho * self._vinfmag**2 + self.pinf
                    self._pressures.append(pressure)
                    self.pressures.append(pressure)
                    self.ports[i]["P"] = pressure
                    
             
                
                

    def _exportPorts(self):
        
        # Specify the CSV file path
        csv_file_path = self.study_directory + '\port_locations.csv'

        # Write data to CSV file
        with open(csv_file_path, 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)

            # Write header
            writer.writerow(['x', 'y', 'z'])

            # Write data rows
            for i in self.portsKeys:
                writer.writerow([self.ports[i]["position"][0], self.ports[i]["position"][1], self.ports[i]["position"][2]])
        csv_file.close()

    def _importPortPressures(self):
        csv_file_path = self.study_directory + "/port_pressures.csv"
        x_points = []
        y_points = []
        z_points = []
        CP = []
        with open(csv_file_path, 'r') as csv_file:
            reader = csv.reader(csv_file)
            # lineCount = len(list(reader))-1
            lineCount = len(self.portsKeys) + 1
            # Read data rows
            line = 0
            for row in reader:
                if line == 0: 
                    line +=1
   
                elif line == lineCount:
                    break
                else:
                    
                    x_points.append(float(row[0]))
                    y_points.append(float(row[1]))
                    z_points.append(float(row[2]))
                    V = float(row[23])
                    CP.append(1 - V**2)
                    line +=1
        csv_file.close()
        
        return x_points, y_points, z_points, CP, 



    def _plotSphere(self):
        # create lists, arrays, and counters
        ci = cj = c = 0
        maxU = 50
        maxV = 50
        u = np.linspace(-np.pi/2, np.pi/2  , maxU)
        v = np.linspace(-np.pi/2, np.pi/2 , maxV)
        points = np.zeros([4,maxU*maxV])
        # get points on sphere
        for i in u:
            cj = 0
            for j in v:
                x , y ,z = self.alphaBetaToCart(i, j)
                points[0,c] = x
                points[1,c] = y
                points[2,c] = z
                cj += 1
                c += 1
            ci += 1
        # get cp for each of the points on the sphere
        for i in range(c):
            point = np.array([points[0,i],points[1,i],points[2,i]])
            mag = np.dot(self._vinf,point)
            t = np.arccos(round(mag,self.tol))
            p = 1-9/4*(np.sin(t))**2
            points[3,i] = p
    
        # create figure
        self.fig2 = plt.figure(1)
        self.ax2 = plt.subplot(projection='3d')
        # set color bar 
        cmap = mpl.cm.winter
        norm = mpl.colors.Normalize(vmax=max(points[3,:]),vmin=min(points[3,:]))
        self.fig2.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),label="C_p")
        # plot points
        self.ax2.scatter(xs=points[0,:],ys=points[1,:],zs=points[2,:],c=points[3,:],s=1,cmap=cmap,norm=norm)
        self.ax2.set_box_aspect((np.ptp(points[0,:]), np.ptp(points[1,:]), np.ptp(points[2,:])))
        self.ax2.set_xlabel("x-axis")
        self.ax2.set_ylabel("y-axis")
        self.ax2.set_zlabel("z-axis")
        
    

    def _plotPorts(self):
        for i in self.portsKeys:
            x = self.ports[i]['position'][0]
            y = self.ports[i]['position'][1]
            z = self.ports[i]['position'][2]
            self.ax2.scatter(x,y,z,color='r')

    def genPorts(self,triad,pitchAngle,rollAngle,yawAngle):
        '''generates an orthogonal triad of port positions on a unit sphere based on 3 euler angles pitch, roll, and yaw in degrees'''
        # inputs
        pitchAngle = np.radians(pitchAngle)
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
            name = "point" + str(triad) + str(i+1)
            self.ports[name]=dict()
            self.ports[name]["position"] = p[:,i]


    def _plotMapping(self):
        # create lists, arrays, and counters
        ci = cj = c = 0
        maxU = 20
        maxV = 20
        u = np.linspace(-np.pi/2, np.pi/2  , maxU)
        v = np.linspace(-np.pi/2, np.pi/2 , maxV)
        points = np.zeros([7,maxU*maxV])
        # get points on sphere
        for i in u:
            cj = 0
            for j in v:
                x , y ,z = self.alphaBetaToCart(i, j)
                points[0,c] = x
                points[1,c] = y
                points[2,c] = z
                points[3,c] = i
                points[4,c] = j
                points[6,c] = self.func([i,j], self.params[0], self.params[1], self.params[2], self.params[3])
                cj += 1
                c += 1
            ci += 1

        # get cp for each of the points on the sphere
        for i in range(c):
            point = np.array([points[0,i],points[1,i],points[2,i]])
            mag = np.dot(self._vinf,point)
            t = np.arccos(round(mag,self.tol))
            cp = 1-9/4*(np.sin(t))**2
            pressure = cp * 0.5 * self.rho * self._vinfmag**2 + self.pinf
            points[5,i] = pressure

        # create figure
        self.fig3 = plt.figure(1)
        self.ax3 = plt.subplot(projection='3d')
        
        # plot points
        self.ax3.scatter(xs=points[3,:],ys=points[4,:],zs=points[5,:],c=points[5,:])
        self.ax3.plot_trisurf(points[3,:],points[4,:],points[6,:],antialiased=True,label="curve fit")
        # self.ax3.set_box_aspect((np.ptp(points[0,:]), np.ptp(points[1,:]), np.ptp(points[2,:])))
        self.ax3.set_xlabel("Alpha (port location)")
        self.ax3.set_ylabel("Beta (port location)")
        self.ax3.set_zlabel("Pressure")
        

        

if __name__ == "__main__":
    study_directory = "studies/PitotProbe"
    sphere1 = Sphere(study_directory,"input_EURP.json",alpha=5,beta=2.5)
    sphere1.run()
