import numpy as np
import matplotlib.pyplot as plt
from panel import SubsonicPanel


if __name__=="__main__":

    # Initialize panel
    panel = SubsonicPanel(1.0, 1.0)
    panel.set_doublet_strength([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    panel.set_source_strength([1.0, 0.0, 0.0])

    # Calculate slice of potentials
    x = np.linspace(-3.0, 3.0, 100)
    z = np.linspace(-3.0, 3.0, 100)
    phi_s_dis = np.zeros((100, 100))
    phi_d_dis = np.zeros((100, 100))
    phi_s_anl = np.zeros((100, 100))
    phi_d_anl = np.zeros((100, 100))
    for i, xi in enumerate(x):
        for j, zj in enumerate(z):

            # Discrete potentials
            phi_s_dis[i,j] = panel.calc_discrete_source_potential([xi, 0.0, zj], 20, 20)
            phi_d_dis[i,j] = panel.calc_discrete_doublet_potential([xi, 0.0, zj], 20, 20)

            # Analytic potentials
            phi_s_anl[i,j] = panel.calc_analytic_source_potential([xi, 0.0, zj])
            phi_d_anl[i,j] = panel.calc_analytic_doublet_potential([xi, 0.0, zj])

    # Plot source potentials
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(6.5, 2))
    ax[0].imshow(phi_s_dis.T)
    ax[1].imshow(phi_s_anl.T)
    ax[2].imshow(abs(phi_s_dis-phi_s_anl).T)
    plt.show()

    # Plot doublet potentials
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(6.5, 2))
    ax[0].imshow(phi_d_dis.T)
    ax[1].imshow(phi_d_anl.T)
    ax[2].imshow(abs(phi_d_dis-phi_d_anl).T)
    plt.show()