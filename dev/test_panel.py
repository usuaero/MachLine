import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from panel import Panel
from singularities import Source, Doublet


if __name__=="__main__":

    # Initialize panel
    verts = np.array([[0.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0]])
    p = Panel(verts)
    print("Panel area: ", p.A)

    # Initialize point source and doublet
    s = Source(np.array([0.0, 0.0, 0.0]), 1.0)
    d = Doublet(np.array([0.0, 0.0, 0.0]), 1.0)

    # Declare array of evaluation positions
    Nx = 1003
    Ny = 1003
    x = np.linspace(-10.0, 10.0, Nx)
    y = np.linspace(-10.0, 10.0, Ny)
    phi_s = np.zeros((Nx, Ny))
    phi_d = np.zeros((Nx, Ny))
    phi_d_up = np.zeros((Nx, Ny))
    phi_d_down = np.zeros((Nx, Ny))

    # Loop through evaluation positions
    for i, xi in enumerate(x):
        for j, yj in enumerate(y):

            # Calculate induced potential
            phi_s[i,j] = s.calc_induced_potential([xi, yj, 0.0])
            phi_d[i,j] = d.calc_induced_potential([xi, yj, 0.0])
            phi_d_up[i,j] = d.calc_induced_potential([xi, yj, 1.0])
            phi_d_down[i,j] = d.calc_induced_potential([xi, yj, -1.0])


    fig, ax = plt.subplots(nrows=2, ncols=2)
    c = cm.get_cmap('viridis', 256)

    # Plot source influence
    ax[0,0].imshow(phi_s, interpolation='nearest', cmap=c, vmin=-1.0, vmax=1.0)
    ax[0,0].set_title('Source Influence z=0')

    # Plot doublet influence
    ax[0,1].imshow(phi_d, interpolation='nearest', cmap=c, vmin=-1.0, vmax=1.0)
    ax[0,1].set_title('Doublet Influence z=0')

    # Plot doublet influence
    ax[1,0].imshow(phi_d_up, interpolation='nearest', cmap=c, vmin=-1.0, vmax=1.0)
    ax[1,0].set_title('Doublet Influence z=1')

    # Plot doublet influence
    ax[1,1].imshow(phi_d_down, interpolation='nearest', cmap=c, vmin=-1.0, vmax=1.0)
    ax[1,1].set_title('Doublet Influence z=-1')

    plt.show()