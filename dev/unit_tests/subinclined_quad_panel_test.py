import numpy as np
import matplotlib.pyplot as plt
from panel import SupersonicSubinclinedPanel


def plot_comparison(x, z, phi_dis, phi_anl, sing_type):
    """Plots a comparison of the predicted potentials."""

    # Initialize figure
    fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(6.5, 2))

    # Get analytic range (since this won't be undefined)
    min_phi = np.min(np.min(phi_anl)).item()
    max_phi = np.max(np.max(phi_anl)).item()

    # Plot potentials
    im1 = ax[1].contourf(x, z, phi_anl.T, levels=50, cmap='hot', vmin=min_phi, vmax=max_phi)
    im0 = ax[0].contourf(x, z, phi_dis.T, levels=im1.levels, cmap='hot', vmin=min_phi, vmax=max_phi)

    # Plot fractional difference
    frac_diff = np.log10(abs((phi_dis-phi_anl)/phi_anl)).T
    im2 = ax[3].contourf(x, z, frac_diff, levels=50, cmap='hot', vmin=-2.0, vmax=0.0)

    # Labels
    ax[0].set_title('$\phi_\{0}$ Discrete'.format(sing_type))
    ax[1].set_title('$\phi_\{0}$ Analytic'.format(sing_type))
    ax[3].set_title('$\log|\Delta_{frac}|$')
    ax[0].set_xlabel("$x$")
    ax[0].set_ylabel("$z$")
    ax[1].set_xlabel("$x$")
    ax[1].set_ylabel("$z$")
    ax[3].set_xlabel("$x$")
    ax[3].set_ylabel("$z$")

    # Add colorbars
    fig.colorbar(im1, cax=ax[2])
    ax[2].set_aspect(20)
    fig.colorbar(im2, cax=ax[4])
    ax[4].set_aspect(20)

    plt.show()


if __name__=="__main__":

    # Initialize panel
    verts = np.array([[0.5, -0.5, -0.5, 0.5],
                      [0.5, 0.5, -0.5, -0.5]])
    panel = SupersonicSubinclinedPanel(verts)
    panel.distribute_points(40, 40)

    # Initialize singularity distributions
    mu_params = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    sigma_params = [0.0, 0.0, 1.0]
    panel.set_doublet_strength(mu_params)
    panel.set_source_strength(sigma_params)

    # Calculate slice of potentials
    Nx = 40
    Nz = 40
    x = np.linspace(0.5, 1.5, Nx)
    z = np.linspace(-0.5, 0.5, Nz)
    phi_s_dis = np.zeros((Nx, Nz))
    phi_d_dis = np.zeros((Nx, Nz))
    phi_s_anl = np.zeros((Nx, Nz))
    phi_d_anl = np.zeros((Nx, Nz))
    for i, xi in enumerate(x):
        for j, zj in enumerate(z):

            # Set point
            #P = [xi, 0.0, zj]
            P = [0.75, xi, zj]

            # Discrete potentials
            phi_s_dis[i,j] = panel.calc_discrete_source_potential(P)
            phi_d_dis[i,j] = panel.calc_discrete_doublet_potential(P)

            # Analytic potentials
            phi_s_anl[i,j] = panel.calc_analytic_source_potential(P)
            phi_d_anl[i,j] = panel.calc_analytic_doublet_potential(P)

    # Plot source potentials
    plot_comparison(x, z, phi_s_dis, phi_s_anl, "sigma")

    # Plot doublet potentials
    plot_comparison(x, z, phi_d_dis, phi_d_anl, "mu")

    ## Plot source potentials subset
    #plot_comparison(x[Nx//2:], z, phi_s_dis[Nx//2:,:], phi_s_anl[Nx//2:,:], "sigma")

    ## Plot doublet potentials subset
    #plot_comparison(x[Nx//2:], z, phi_d_dis[Nx//2:,:], phi_d_anl[Nx//2:,:], "mu")