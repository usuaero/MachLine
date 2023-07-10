import numpy as np

if __name__=="__main__":

    # Limits
    x_max = 5.0
    x_min = -1.0
    y_max = 0.0
    y_min = 0.0
    z_max = 1.0
    z_min = -1.0

    # Number of points to distribute
    Nx = 20
    Ny = 1
    Nz = 20

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
    filename = "test/input_files/root_xz_sample_points.csv"
    with open(filename, 'w') as points_file:

        # Header
        print("x,y,z", file=points_file)

        # loop through points
        for i in range(len(X)):
            print("{0},{1},{2}".format(X[i], Y[i], Z[i]), file=points_file)