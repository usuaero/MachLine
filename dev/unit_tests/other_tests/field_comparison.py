import numpy as np
import matplotlib.pyplot as plt
from panel import Panel
from singularities import Source, Doublet


def compare_on_box(verts, lims, N, filename):
    # Compares the panel and point calculations over the specified box

    # Initialize panel
    p = Panel(verts, 1.0, 1.0)

    # Initialize equivalent point singularities
    s = Source(np.array([p.c[0], p.c[1], 0.0]), p.A)
    d = Doublet(np.array([p.c[0], p.c[1], 0.0]), p.A)

    # Declare grid coordinate vectors
    x = np.linspace(lims[0][0], lims[0][1], N[0])
    y = np.linspace(lims[1][0], lims[1][1], N[1])
    z = np.linspace(lims[2][0], lims[2][1], N[2])

    # Get coordinate matrices
    X, Y, Z = np.meshgrid(x, y, z)
    X = X.flatten()
    Y = Y.flatten()
    Z = Z.flatten()

    # Initialize storage
    N = X.size
    phi_s_point = np.zeros(N)
    phi_d_point = np.zeros(N)
    phi_s_panel = np.zeros(N)
    phi_d_panel = np.zeros(N)

    # Loop through points
    for i in range(N):

        # Evaluate
        P = np.array([X[i], Y[i], Z[i]])
        phi_s_point[i] = s.calc_induced_potential(P)
        phi_d_point[i] = d.calc_induced_potential(P)
        phi_s_panel[i] = p.calc_induced_source_potential(P)
        phi_d_panel[i] = p.calc_induced_doublet_potential(P)

    # Write to file
    with open(filename, 'w') as file_handle:

        # Header
        print('x,y,z,phi_s_point,phi_d_point,phi_s_panel,phi_d_panel', file=file_handle)

        # Lines
        for i in range(N):

            line = ','.join(map(str, [X[i], Y[i], Z[i], phi_s_point[i], phi_d_point[i], phi_s_panel[i], phi_d_panel[i]]))
            print(line, file=file_handle)


if __name__=="__main__":

    # Box comparison
    verts = np.array([[0.0, 1.0, 2.0, 1.0],
                      [0.0, -0.5, 0.0, 0.5]])
    compare_on_box(verts, [[-1.0, 10.0], [-10.0, 10.0], [-10.0, 10.0]], [50, 100, 100], "dev/unit_tests/test.csv")