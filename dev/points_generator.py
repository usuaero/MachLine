import numpy as np

if __name__=="__main__":

    # Limits
    x_max = 1.05
    x_min = 0.95
    y_max = 0.3
    y_min = -0.00
    z_max = 0.05
    z_min = -0.05

    # Number of points to distribute
    Nx = 30
    Ny = 90
    Nz = 30

    # Distributions in each direction
    x = np.linspace(x_min, x_max, Nx)
    y = np.linspace(y_min, y_max, Ny)
    z = np.linspace(z_min, z_max, Nz)

    # Grid
    X, Y, Z = np.meshgrid(x, y, z)
    X = X.flatten()
    Y = Y.flatten()
    Z = Z.flatten()

    # write to file
    filename = "dev/delta_wing_trailing_edge_sample_points.csv"
    with open(filename, 'w') as points_file:

        # Header
        print("x,y,z", file=points_file)

        # loop through points
        for i in range(len(X)):
            print("{0},{1},{2}".format(X[i], Y[i], Z[i]), file=points_file)