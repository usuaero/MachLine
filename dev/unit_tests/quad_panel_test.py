import numpy as np
import matplotlib.pyplot as plt
from panel import SubsonicPanel
from mpl_toolkits.axes_grid1 import make_axes_locatable


if __name__=="__main__":

    # Initialize panel
    panel = SubsonicPanel(1.0, 1.0)
    mu_params = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    sigma_params = [0.0, 0.0, 1.0]
    panel.set_doublet_strength(mu_params)
    panel.set_source_strength(sigma_params)

    # Calculate slice of potentials
    x = np.linspace(-1.0, 1.0, 30)
    z = np.linspace(-1.0, 1.0, 30)
    phi_s_dis = np.zeros((30, 30))
    phi_d_dis = np.zeros((30, 30))
    phi_s_anl = np.zeros((30, 30))
    phi_d_anl = np.zeros((30, 30))
    for i, xi in enumerate(x):
        for j, zj in enumerate(z):

            # Set point
            P = [xi, xi, zj]

            # Discrete potentials
            phi_s_dis[i,j] = panel.calc_discrete_source_potential(P, 50, 50)
            phi_d_dis[i,j] = panel.calc_discrete_doublet_potential(P, 50, 50)

            # Analytic potentials
            phi_s_anl[i,j] = panel.calc_analytic_source_potential(P)
            phi_d_anl[i,j] = panel.calc_analytic_doublet_potential(P)

    # Plot source potentials
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(6.5, 2))
    min_phi = min(np.min(np.min(phi_s_dis)).item(), np.min(np.min(phi_s_anl)).item())
    max_phi = max(np.max(np.max(phi_s_dis)).item(), np.max(np.max(phi_s_anl)).item())
    im0 = ax[0].imshow(phi_s_dis.T, cmap='hot', vmin=min_phi, vmax=max_phi)
    im1 = ax[1].imshow(phi_s_anl.T, cmap='hot', vmin=min_phi, vmax=max_phi)
    im2 = ax[2].imshow(np.log10(abs((phi_s_dis-phi_s_anl)/phi_s_anl)).T, cmap='hot')
    ax[0].set_title('Discrete')
    ax[1].set_title('Analytic')
    ax[2].set_title('$\log|\Delta_{frac}|$')
    val_cbar_ax = fig.add_axes([0.62, 0.15, 0.01, 0.6])
    fig.colorbar(im0, cax=val_cbar_ax)
    err_cbar_ax = fig.add_axes([0.93, 0.15, 0.01, 0.6])
    fig.colorbar(im2, cax=err_cbar_ax)
    fig.suptitle('Source-Induced Potential for {0}'.format(sigma_params))
    plt.show()

    # Plot doublet potentials
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(6.5, 2))
    min_phi = min(np.min(np.min(phi_d_dis)).item(), np.min(np.min(phi_d_anl)).item())
    max_phi = max(np.max(np.max(phi_d_dis)).item(), np.max(np.max(phi_d_anl)).item())
    im0 = ax[0].imshow(phi_d_dis.T, cmap='hot', vmin=min_phi, vmax=max_phi)
    im1 = ax[1].imshow(phi_d_anl.T, cmap='hot', vmin=min_phi, vmax=max_phi)
    im2 = ax[2].imshow(np.log10(abs((phi_d_dis-phi_d_anl)/phi_d_anl)).T, cmap='hot')
    ax[0].set_title('Discrete')
    ax[1].set_title('Analytic')
    ax[2].set_title('$\log|\Delta_{frac}|$')
    val_cbar_ax = fig.add_axes([0.62, 0.15, 0.01, 0.6])
    fig.colorbar(im0, cax=val_cbar_ax)
    err_cbar_ax = fig.add_axes([0.93, 0.15, 0.01, 0.6])
    fig.colorbar(im2, cax=err_cbar_ax)
    fig.suptitle('Doublet-Induced Potential for {0}'.format(mu_params))
    plt.show()