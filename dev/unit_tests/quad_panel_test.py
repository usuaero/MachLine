import numpy as np
import matplotlib.pyplot as plt
from panel import SubsonicPanel
from mpl_toolkits.axes_grid1 import make_axes_locatable


if __name__=="__main__":

    # Initialize panel
    panel = SubsonicPanel(1.0, 1.0)
    panel.set_doublet_strength([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    panel.set_source_strength([1.0, 0.0, 0.0])

    # Calculate slice of potentials
    x = np.linspace(-1.0, 1.0, 30)
    z = np.linspace(-1.0, 1.0, 30)
    phi_s_dis = np.zeros((30, 30))
    phi_d_dis = np.zeros((30, 30))
    phi_s_anl = np.zeros((30, 30))
    phi_d_anl = np.zeros((30, 30))
    for i, xi in enumerate(x):
        for j, zj in enumerate(z):

            # Discrete potentials
            phi_s_dis[i,j] = panel.calc_discrete_source_potential([xi, 0.0, zj], 40, 40)
            phi_d_dis[i,j] = panel.calc_discrete_doublet_potential([xi, 0.0, zj], 40, 40)

            # Analytic potentials
            phi_s_anl[i,j] = panel.calc_analytic_source_potential([xi, 0.0, zj])
            phi_d_anl[i,j] = panel.calc_analytic_doublet_potential([xi, 0.0, zj])

    # Plot source potentials
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(6.5, 2))
    im0 = ax[0].imshow(phi_s_dis.T, cmap='hot')
    im1 = ax[1].imshow(phi_s_anl.T, cmap='hot')
    im2 = ax[2].imshow(np.log10(abs(phi_s_dis-phi_s_anl)).T, cmap='hot')
    ax[0].set_title('Discrete')
    ax[1].set_title('Analytic')
    ax[2].set_title('$\log|\Delta|$')
    fig.colorbar(im2)
    fig.suptitle('Source-Induced Potential')
    plt.show()

    # Plot doublet potentials
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(6.5, 2))
    im0 = ax[0].imshow(phi_d_dis.T, cmap='hot')
    im1 = ax[1].imshow(phi_d_anl.T, cmap='hot')
    im2 = ax[2].imshow(np.log10(abs(phi_d_dis-phi_d_anl)).T, cmap='hot')
    ax[0].set_title('Discrete')
    ax[1].set_title('Analytic')
    ax[2].set_title('$\log|\Delta|$')
    fig.colorbar(im2)
    fig.suptitle('Doublet-Induced Potential')
    plt.show()