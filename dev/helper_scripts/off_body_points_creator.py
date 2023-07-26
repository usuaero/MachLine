import numpy as np
import csv

def point_block(x_min=-0.25,x_max=1.5,y_min=-1.8,y_max=-1.8,z_min=-0.25,z_max=0.25,Nx=200,Ny=1,Nz=100):
    # Create point block
    #x_min = -1.0
    #x_max = 3.0
    #y_min = -2.0
    #y_max = 2.0
    #z_min = 0.001
    #z_max = 0.001

    points = []

    x = np.linspace(x_min,x_max,Nx)
    y = np.linspace(y_min,y_max,Ny)
    z = np.linspace(z_min,z_max,Nz)
    
    X,Y,Z = np.meshgrid(x,y,z)

    for i in range(len(X)):
        for j in range(len(X[i])):
            for k in range(len(X[i][j])):
                points.append([X[i][j][k],Y[i][j][k],Z[i][j][k]])

    return points


def sphere_surface_slice(center,radius,offset=1e-6,n_points=100):
    thetas = np.linspace(0,np.pi,n_points)

    points = []

    for i, theta in enumerate(thetas):
        x = (radius+offset)*np.cos(theta) + center[0]
        y = (radius+offset)*np.sin(theta) + center[1]
        z = center[2]

        points.append([x,y,z])

    return points, thetas


def cone_ray_slice(tip, ray_angle, n_points=100):
    x_s = np.linspace(tip[0],0,n_points)

    points = []

    for i, x in enumerate(x_s):
        if i == 0:
            pass
        else:
            y = (tip[0]-x)*np.tan(np.deg2rad(ray_angle)) + tip[1]
            z = tip[2]

            points.append([x,y,z])
    
    return points

def write_points(filename,points):
    
    header = ['x','y','z']

    with open(filename, 'w') as csvfile:

        csvwriter = csv.writer(csvfile)

        csvwriter.writerow(header)

        csvwriter.writerows(points)


if __name__ =="__main__":
    points = point_block()

    write_points("studies/input_files/diamond_5_deg_sample_points.csv",points)

