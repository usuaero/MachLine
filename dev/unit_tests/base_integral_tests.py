# Checks the calculation of hH(1,1,3) and F(1,1,1) against Gaussian quadrature for subinclined, supersonic panels

import numpy as np
import scipy.integrate as sint
from dev.unit_tests.panel import SupersonicSubinclinedPanel

P = np.array([5.0, 1.0, 1.0])


def hH113_integrand(y, x):
    # Returns the integrand of hH(1,1,1) for P

    R = np.sqrt((x-P[0])**2 - (y-P[1])**2 - P[2]**2)
    return P[2]/R**3


def F111_integrand(l):
    # Returns the integrand of F(1,1,1) for P

    # First edge
    #R = np.sqrt((l-P[0])**2 - (1.0-P[1])**2 - P[2]**2)

    # Second edge
    #R = np.sqrt((-1.0-P[0])**2 - (l-P[1])**2 - P[2]**2)

    # Third edge
    #R = np.sqrt((l-P[0])**2 - (-1.0-P[1])**2 - P[2]**2)

    # Fourth edge
    R = np.sqrt((1.0-P[0])**2 - (l-P[1])**2 - P[2]**2)
    print(R)

    return 1.0/R


if __name__=="__main__":

    # Initialize panel
    verts = np.array([[1.0, -1.0, -1.0, 1.0],
                      [1.0, 1.0, -1.0, -1.0]])
    panel = SupersonicSubinclinedPanel(verts)

    # Get integrals
    geom = panel.calc_geom(P)
    I = panel.calc_F_integrals(geom)
    panel.calc_H_integrals(geom, I)

    # Calculate F(1,1,1) using quadrature
    F111_quad = sint.quad(F111_integrand, -1.0, 1.0)

    # Calculate hH(1,1,3) using quadrature
    hH113_quad = sint.dblquad(hH113_integrand, -1.0, 1.0, -1.0, 1.0)

    # Compare
    print()
    print("F(1,1,1)")
    print("    Analytic:", I.F111)
    print("    Quadrature: ", F111_quad)
    print()
    print("hH(1,1,3)")
    print("    Analytic:", I.hH113)
    print("    Quadrature: ", hH113_quad)