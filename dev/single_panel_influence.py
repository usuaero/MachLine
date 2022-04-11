import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from panel import Panel
from singularities import Source, Doublet


if __name__=="__main__":

    # Initialize panel
    verts = np.array([[0.0, 1.0, 2.0, 1.0],
                      [0.0, -0.5, 0.0, 0.5]])
    p = Panel(verts, 1.0, 1.0)

    # Initialize equivalent point singularities
    s = Source(np.array([p.c[0], p.c[1], 0.0]), p.A)
    d = Doublet(np.array([p.c[0], p.c[1], 0.0]), p.A)

    # Initialize ray of evaluation points
    P0 = np.array([4.0, 0.0, -1.0])
    d0 = np.array([1.0, 0.0, 0.0])
    d0 /= np.linalg.norm(d0)
    k = np.linspace(0.0, 10.0, 500)
    phi_s_point = np.zeros_like(k)
    phi_d_point = np.zeros_like(k)
    phi_s_panel = np.zeros_like(k)
    phi_d_panel = np.zeros_like(k)

    # Calculate induced potentials
    for i, ki in enumerate(k):
        
        # Get point
        P = P0 + d0*ki

        # Evaluate
        phi_s_point[i] = s.calc_induced_potential(P)
        phi_d_point[i] = d.calc_induced_potential(P)
        phi_s_panel[i] = p.calc_induced_source_potential(P)
        phi_d_panel[i] = p.calc_induced_doublet_potential(P)

    # Plot
    fig, ax = plt.subplots(ncols=2, figsize=(12.0, 7.0))
    plt.suptitle("Panel with purely subsonic edges.")

    # Source influences
    ax[0].plot(k, phi_s_point, 'k--', label='Point')
    ax[0].plot(k, phi_s_panel, 'k-', label='Panel')
    ax[0].set_xlabel('Distance')
    ax[0].set_ylabel('$\phi_s$')
    ax[0].legend()

    # Doublet influences
    ax[1].plot(k, phi_d_point, 'k--', label='Point')
    ax[1].plot(k, phi_d_panel, 'k-', label='Panel')
    ax[1].set_xlabel('Distance')
    ax[1].set_ylabel('$\phi_d$')
    ax[1].legend()

    plt.show()