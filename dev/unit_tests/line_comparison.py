import numpy as np
import matplotlib.pyplot as plt
from panel import Panel
from singularities import Source, Doublet


def plot_comparison(verts, P0, d0, figname=None):

    # Initialize panel
    p = Panel(verts, 1.0, 1.0)

    # Initialize equivalent point singularities
    s = Source(np.array([p.c[0], p.c[1], 0.0]), p.A)
    d = Doublet(np.array([p.c[0], p.c[1], 0.0]), p.A)

    # Initialize ray of evaluation points
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

    if figname is not None:
        plt.savefig(figname)
        plt.close()
    else:
        plt.show()


if __name__=="__main__":

    # Panel with subsonic edges
    verts = np.array([[0.0, 1.0, 2.0, 1.0],
                      [0.0, -0.5, 0.0, 0.5]])

    # Run line in z direction
    P0 = np.array([10.0, 0.0, -5.0])
    d0 = np.array([0.0, 0.0, 1.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/subsonic_edges_z_trans_totally_in.pdf')

    # Run lin in x direction
    P0 = np.array([3.5, 0.0, 1.0])
    d0 = np.array([1.0, 0.0, 0.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/subsonic_edges_x_trans_totally_in.pdf')

    # Panel with supersonic edges
    verts = np.array([[0.0, 1.0, 2.0, 1.0],
                      [0.0, -1.5, 0.0, 1.5]])
    
    # Run line in z direction
    P0 = np.array([10.0, 0.0, -5.0])
    d0 = np.array([0.0, 0.0, 1.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/supersonic_edges_z_trans_totally_in.pdf')

    # Run line in x direction
    P0 = np.array([3.5, 0.0, 1.0])
    d0 = np.array([1.0, 0.0, 0.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/supersonic_edges_x_trans_totally_in.pdf')

    # Triangular panel
    verts = np.array([[0.0, 0.0, 1.0],
                      [1.0, -1.0, 0.0]])

    # Run line through Mach wedge region
    P0 = np.array([0.9, 0.0, -5.0])
    d0 = np.array([0.0, 0.0, 1.0])
    plot_comparison(verts, P0, d0)