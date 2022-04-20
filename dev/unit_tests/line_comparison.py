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
    N_line = 500
    d0 /= np.linalg.norm(d0)
    k = np.linspace(0.0, 10.0, N_line)

    # Initialize storage
    phi_s_point = np.zeros_like(k)
    phi_d_point = np.zeros_like(k)
    phi_s_panel = np.zeros_like(k)
    phi_d_panel = np.zeros_like(k)
    hH113 = np.zeros_like(k)
    F111 = np.zeros((p.N, N_line))
    a = np.zeros((p.N, N_line))

    # Calculate induced potentials
    for i, ki in enumerate(k):
        
        # Get point
        P = P0 + d0*ki

        # Evaluate
        phi_s_point[i] = s.calc_induced_potential(P)
        phi_d_point[i] = d.calc_induced_potential(P)
        phi_s_panel[i], hH113[i], F111[:,i], a[:,i] = p.calc_induced_source_potential(P)
        phi_d_panel[i] = p.calc_induced_doublet_potential(P)

    # Plot
    fig, ax = plt.subplots(ncols=2, figsize=(12.0, 7.0))

    # Source influences
    ax[0].plot(k, phi_s_point, 'k--', label='Point')
    ax[0].plot(k, phi_s_panel, 'k-', label='Panel')
    ax[0].set_xlabel('Distance')
    ax[0].set_ylabel('$\phi_s$')
    ax[0].legend()
    if np.min(phi_s_panel) > 0.0:
        ax[0].set_ylim(bottom=0.0, top=np.max(phi_s_panel)*1.2)
    elif np.min(phi_s_panel) == 0.0:
        ax[0].set_ylim(bottom=-0.1, top=np.max(phi_s_panel)*1.2)
    else:
        ax[0].set_ylim(bottom=np.min(phi_s_panel)*1.2, top=np.max(phi_s_panel)*1.2)

    # Doublet influences
    ax[1].plot(k, phi_d_point, 'k--', label='Point')
    ax[1].plot(k, phi_d_panel, 'k-', label='Panel')
    ax[1].set_xlabel('Distance')
    ax[1].set_ylabel('$\phi_d$')
    ax[1].legend()
    ax[1].set_ylim(bottom=np.min(phi_d_panel)*1.2, top=np.max(phi_d_panel)*1.2)

    if figname is not None:
        plt.savefig(figname)
        plt.close()
    else:
        plt.show()

    # Plot
    fig, ax = plt.subplots(ncols=3, figsize=(12.0, 5.0))

    # hH(1,1,3)
    ax[0].plot(k, hH113, 'k-')
    ax[0].set_xlabel('Distance')
    ax[0].set_ylabel('$hH(1,1,3)$')

    # F111
    for i in range(p.N):
        ax[1].plot(k, F111[i,:], label=str(i))
    ax[1].set_xlabel('Distance')
    ax[1].set_ylabel('$F(1,1,1)$')
    ax[1].legend()

    # a
    for i in range(p.N):
        ax[2].plot(k, a[i,:], label=str(i))
    ax[2].set_xlabel('Distance')
    ax[2].set_ylabel('$a$')
    ax[2].legend()

    if figname is not None:
        plt.savefig(figname.replace('.pdf', '_integrals.pdf'))
        plt.close()
    else:
        plt.show()


if __name__=="__main__":

    # Triangular panel
    verts = np.array([[0.0, 0.0, 1.0],
                      [1.0, -1.0, 0.0]])

    # Run line through Mach wedge region
    P0 = np.array([0.9, 0.0, -5.0])
    d0 = np.array([0.0, 0.0, 1.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/mach_wedge.pdf')



    # Panel with subsonic edges
    verts = np.array([[0.0, 1.0, 2.0, 1.0],
                      [0.0, -0.5, 0.0, 0.5]])

    # Run line in z direction
    P0 = np.array([10.0, 0.0, -5.0])
    d0 = np.array([0.0, 0.0, 1.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/subsonic_edges_z_trans_totally_in.pdf')
    
    # Run line in z direction passing over Mach cone
    P0 = np.array([3, 0.0, -5.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/subsonic_edges_z_trans_intersecting.pdf')

    # Run line in x direction
    P0 = np.array([3, 0.0, 1.0])
    d0 = np.array([1.0, 0.0, 0.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/subsonic_edges_x_trans_totally_in.pdf')
    
    # Run line in y direction passing over Mach cone
    P0 = np.array([3, -5.0, 1.0])
    d0 = np.array([0.0, 1.0, 0.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/subsonic_edges_y_trans_intersecting.pdf')



    # Panel with supersonic edges
    verts = np.array([[0.0, 1.0, 2.0, 1.0],
                      [0.0, -1.5, 0.0, 1.5]])
    
    # Run line in z direction
    P0 = np.array([10.0, 0.0, -5.0])
    d0 = np.array([0.0, 0.0, 1.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/supersonic_edges_z_trans_totally_in.pdf')
    
    # Run line in z direction passing over Mach cone
    P0 = np.array([3, 0.0, -5.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/supersonic_edges_z_trans_intersecting.pdf')

    # Run line in x direction
    P0 = np.array([3, 0.0, 1.0])
    d0 = np.array([1.0, 0.0, 0.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/supersonic_edges_x_trans_totally_in.pdf')
    
    # Run line in y direction passing over Mach cone
    P0 = np.array([3, -5.0, 1.0])
    d0 = np.array([0.0, 1.0, 0.0])
    plot_comparison(verts, P0, d0, 'dev/unit_tests/supersonic_edges_y_trans_intersecting.pdf')