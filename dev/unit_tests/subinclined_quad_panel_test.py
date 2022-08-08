import numpy as np
import matplotlib.pyplot as plt
from panel import SupersonicSubinclinedPanel


if __name__=="__main__":

    # Initialize panel
    verts = np.array([[0.5, -1.5, -0.5, 0.5],
                      [0.5, 0.5, -0.5, -0.5]])
    panel = SupersonicSubinclinedPanel(verts)
    panel.distribute_points(40, 40)

    # Initialize singularity distributions
    mu_params = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    sigma_params = [0.0, 1.0, 0.0]
    panel.set_doublet_strength(mu_params)
    panel.set_source_strength(sigma_params)

    # Calculate slice of potentials
    Nx = 60
    Nz = 30
    x = np.linspace(0.0, 4.0, Nx)
    z = np.linspace(-1.0, 1.0, Nz)
    phi_s_dis = np.zeros((Nx, Nz))
    phi_d_dis = np.zeros((Nx, Nz))
    phi_s_anl = np.zeros((Nx, Nz))
    phi_d_anl = np.zeros((Nx, Nz))
    for i, xi in enumerate(x):
        for j, zj in enumerate(z):

            # Set point
            #P = [xi, 0.0, zj]
            P = [xi, xi-1.5, zj]

            # Discrete potentials
            phi_s_dis[i,j] = panel.calc_discrete_source_potential(P)
            phi_d_dis[i,j] = panel.calc_discrete_doublet_potential(P)

            # Analytic potentials
            phi_s_anl[i,j] = panel.calc_analytic_source_potential(P)
            phi_d_anl[i,j] = panel.calc_analytic_doublet_potential(P)

    # Plot source potentials
    fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(6.5, 2))

    min_phi = np.min(np.min(phi_s_anl)).item()
    max_phi = np.max(np.max(phi_s_anl)).item()
    im0 = ax[0].imshow(phi_s_dis.T, cmap='hot', vmin=min_phi, vmax=max_phi)
    im1 = ax[1].imshow(phi_s_anl.T, cmap='hot', vmin=min_phi, vmax=max_phi)

    frac_diff = np.log10(abs((phi_s_dis-phi_s_anl)/phi_s_anl)).T
    im2 = ax[3].imshow(frac_diff, cmap='hot', vmin=-2.0, vmax=0.0)

    ax[0].set_title('$\phi_\sigma$ Discrete')
    ax[1].set_title('$\phi_\sigma$ Analytic')
    ax[3].set_title('$\log|\Delta_{frac}|$')

    fig.colorbar(im0, cax=ax[2])
    ax[2].set_aspect(10)
    fig.colorbar(im2, cax=ax[4])
    ax[4].set_aspect(10)

    plt.show()

    # Plot doublet potentials
    fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(6.5, 2))

    min_phi = np.min(np.min(phi_d_anl)).item()
    max_phi = np.max(np.max(phi_d_anl)).item()
    im0 = ax[0].imshow(phi_d_dis.T, cmap='hot', vmin=min_phi, vmax=max_phi)
    im1 = ax[1].imshow(phi_d_anl.T, cmap='hot', vmin=min_phi, vmax=max_phi)

    frac_diff = np.log10(abs((phi_d_dis-phi_d_anl)/phi_d_anl)).T
    im2 = ax[3].imshow(frac_diff, cmap='hot', vmin=-2.0, vmax=0.0)

    ax[0].set_title('$\phi_\mu$ Discrete')
    ax[1].set_title('$\phi_\mu$ Analytic')
    ax[3].set_title('$\log|\Delta_{frac}|$')

    fig.colorbar(im0, cax=ax[2])
    ax[2].set_aspect(10)
    fig.colorbar(im2, cax=ax[4])
    ax[4].set_aspect(10)

    plt.show()

    # Plot source potentials subset
    fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(6.5, 2))

    min_phi = np.min(np.min(phi_s_anl[Nx//2:,:])).item()
    max_phi = np.max(np.max(phi_s_anl[Nx//2:,:])).item()
    im0 = ax[0].imshow(phi_s_dis[Nx//2:,:].T, cmap='hot', vmin=min_phi, vmax=max_phi)
    im1 = ax[1].imshow(phi_s_anl[Nx//2:,:].T, cmap='hot', vmin=min_phi, vmax=max_phi)

    frac_diff = np.log10(abs((phi_s_dis[Nx//2:,:]-phi_s_anl[Nx//2:,:])/phi_s_anl[Nx//2:,:])).T
    im2 = ax[3].imshow(frac_diff, cmap='hot', vmin=-2.0, vmax=0.0)

    ax[0].set_title('$\phi_\sigma$ Discrete')
    ax[1].set_title('$\phi_\sigma$ Analytic')
    ax[3].set_title('$\log|\Delta_{frac}|$')

    fig.colorbar(im0, cax=ax[2])
    ax[2].set_aspect(10)
    fig.colorbar(im2, cax=ax[4])
    ax[4].set_aspect(10)

    plt.show()

    # Plot doublet potentials subset
    fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(6.5, 2))

    min_phi = np.min(np.min(phi_d_anl[Nx//2:,:])).item()
    max_phi = np.max(np.max(phi_d_anl[Nx//2:,:])).item()
    im0 = ax[0].imshow(phi_d_dis[Nx//2:,:].T, cmap='hot', vmin=min_phi, vmax=max_phi)
    im1 = ax[1].imshow(phi_d_anl[Nx//2:,:].T, cmap='hot', vmin=min_phi, vmax=max_phi)

    frac_diff = np.log10(abs((phi_d_dis[Nx//2:,:]-phi_d_anl[Nx//2:,:])/phi_d_anl[Nx//2:,:])).T
    im2 = ax[3].imshow(frac_diff, cmap='hot', vmin=-2.0, vmax=0.0)

    ax[0].set_title('$\phi_\mu$ Discrete')
    ax[1].set_title('$\phi_\mu$ Analytic')
    ax[3].set_title('$\log|\Delta_{frac}|$')

    fig.colorbar(im0, cax=ax[2])
    ax[2].set_aspect(10)
    fig.colorbar(im2, cax=ax[4])
    ax[4].set_aspect(10)

    plt.show()