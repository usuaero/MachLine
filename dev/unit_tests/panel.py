import numpy as np


def inner2(x,y):
    return x[0]*y[0] + x[1]*y[1]

K_inv = 1.0/(4.0*np.pi)


class Geometry:
    """A class containing all geometric parameters for a point-panel pair."""

    def __init__(self):
        self.a = np.zeros(4)
        self.g2 = np.zeros(4)
        self.l1 = np.zeros(4)
        self.l2 = np.zeros(4)
        self.R1 = np.zeros(4)
        self.R2 = np.zeros(4)
        self.v_xi = np.zeros(4)
        self.v_eta = np.zeros(4)
        self.h = 0.0
        self.h2 = 0.0
        self.b = np.ones(4)


class Integrals:
    """A class containing all necessary ints."""

    def __init__(self):
        self.hH113 = 0.0
        self.H111 = 0.0
        self.H211 = 0.0
        self.H121 = 0.0
        self.H213 = 0.0
        self.H123 = 0.0
        self.H313 = 0.0
        self.H133 = 0.0
        self.H223 = 0.0
        self.F111 = np.zeros(4)
        self.F211 = np.zeros(4)
        self.F121 = np.zeros(4)


class SubsonicPanel:
    """A class defining a rectangular, subsonic, quadratic-doublet-linear-source panel in incompressible flow.
    The panel lies in the x-y (z=0) plane and is centered at the origin.

           0
       ----------
       |        |
     1 |        | 3
       |        |
       ----------
            2
    
    Parameters
    ----------
    x_dim : float
        Width in the x direction.

    y_dim : float
        Width in the y direction.
    """

    def __init__(self, x_dim, y_dim):

        # Store
        self.x_dim = x_dim
        self.y_dim = y_dim
        self.verts = np.zeros((2,4))
        
        # First vertex
        self.verts[0,0] = 0.5*x_dim
        self.verts[1,0] = 0.5*y_dim
        
        # Second vertex
        self.verts[0,1] = -0.5*x_dim
        self.verts[1,1] = 0.5*y_dim
        
        # Third vertex
        self.verts[0,2] = -0.5*x_dim
        self.verts[1,2] = -0.5*y_dim
        
        # Fourth vertex
        self.verts[0,3] = 0.5*x_dim
        self.verts[1,3] = -0.5*y_dim

        # Initialize a few things
        self.mu_params = np.zeros(6)
        self.sigma_params = np.zeros(3)


    def set_doublet_strength(self, mu_params):
        self.mu_params = mu_params


    def set_source_strength(self, sigma_params):
        self.sigma_params = sigma_params


    def calc_geom(self, P):
        """Calculates the geometry for the panel-point pair.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        geom : Geometry
            Geometry of the point relative to the panel
        """

        # Initialize
        geom = Geometry()
        geom.h = P[2]
        geom.h2 = geom.h**2

        # Outward normals
        geom.v_eta[0] = 1.0
        geom.v_eta[2] = -1.0
        geom.v_xi[1] = 1.0
        geom.v_xi[3] = -1.0

        # Loop through edges
        for i in range(4):

            # Displacement to first vertex
            d = self.verts[:,i] - P[0:2]
            
            # Perpendicular in-plane distance
            geom.a[i] = d[0]*geom.v_xi[i] + d[1]*geom.v_eta[i]

            # Tangential in-plane distance
            geom.l1[i] = -d[0]*geom.v_eta[i] + d[1]*geom.v_xi[i]

            # Displacement to second vertex
            d = self.verts[:,(i+1)%4] - P[0:2]

            # Tangential in-plane distance
            geom.l2[i] = -d[0]*geom.v_eta[i] + d[1]*geom.v_xi[i]

        # Perpendicular distances
        geom.g2 = geom.a**2 + geom.h2

        # Radial distances
        geom.R1 = np.sqrt(geom.l1**2 + geom.g2)
        geom.R2 = np.sqrt(geom.l2**2 + geom.g2)

        return geom


    def calc_F_integrals(self, geom):
        """Calculates necessary F integrals.
        
        Parameters
        ----------
        geom : Geometry
            Geometry of the point relative to the panel

        Returns
        -------
        ints : Integrals
            Container of F ints.
        """

        ints = Integrals()

        # Loop through edges
        for i in range(4):

            # Within edge
            if geom.l1[i]*geom.l2[i] < 0.0:

                ints.F111[i] = np.log(((geom.R1[i] - geom.l1[i])*(geom.R2[i] + geom.l2[i]))/geom.g2[i])

            # Outside edge
            else:

                ints.F111[i] = np.sign(geom.l1[i])*np.log((geom.R2[i]+np.abs(geom.l2[i]))/(geom.R1[i]+np.abs(geom.l1[i])))

        # Calculate F(2,1,1) and F(1,2,1)
        # These are correct because the linear source calculations are correct
        ints.F121 = geom.a*geom.v_eta*ints.F111 + geom.v_xi*(geom.R2-geom.R1)
        ints.F211 = geom.a*geom.v_xi*ints.F111 - geom.v_eta*(geom.R2-geom.R1)

        return ints


    def calc_H_integrals(self, geom, ints):
        """Calculates H(M,N,K) integrals.
        
        Parameters
        ----------
        geom : Geometry
            Geometry of the point relative to the panel

        ints : Integrals
            Integrals calculated thus far
        """

        # Loop through edges
        ints.hH113 = 0.0
        for i in range(4):
            
            # Intermediate quantities
            c1 = geom.g2[i] + np.abs(geom.h)*geom.R1[i]
            c2 = geom.g2[i] + np.abs(geom.h)*geom.R2[i]

            # Integral for edge
            S = geom.a[i]*(geom.l2[i]*c1 - geom.l1[i]*c2)
            C = c1*c2 + geom.a[i]**2*geom.l1[i]*geom.l2[i]

            # Sum
            ints.hH113 += np.arctan2(S, C)

        # Apply sign factor
        ints.hH113 *= np.sign(geom.h)

        # Calculate H(1,1,1)
        ints.H111 = -geom.h*ints.hH113 + np.sum(geom.a*ints.F111).item()

        # Calcualte H(2,1,3) and H(1,2,3)
        ints.H213 = -np.sum(geom.v_xi*ints.F111).item()
        ints.H123 = -np.sum(geom.v_eta*ints.F111).item()
        
        # Calculate H(2,1,1) and H(1,2,1)
        ints.H211 = 0.5*(-geom.h2*ints.H213 + np.sum(geom.a*ints.F211).item())
        ints.H121 = 0.5*(-geom.h2*ints.H123 + np.sum(geom.a*ints.F121).item())

        # Calculate H(3,1,3), H(2,2,3), and H(1,3,3)
        ints.H313 = np.sum(geom.v_eta*ints.F121).item() - geom.h*ints.hH113
        ints.H223 = -np.sum(geom.v_xi*ints.F121).item()
        ints.H133 = ints.H111 - np.sum(geom.v_eta*ints.F121).item()


    def calc_analytic_source_potential(self, P):
        """Calculates the potential induced assuming a continuous distribution of source strength.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        phi_s : float
            source-induced potential.
        """

        # Calculate geometry
        P = np.array(P)
        geom = self.calc_geom(P)

        # Calculate necessary ints
        I = self.calc_F_integrals(geom)
        self.calc_H_integrals(geom, I)

        phi_s = self.sigma_params[0]*I.H111 + self.sigma_params[1]*(I.H111*P[0] + I.H211) + self.sigma_params[2]*(I.H111*P[1] + I.H121)
        return -phi_s*K_inv


    def calc_analytic_doublet_potential(self, P):
        """Calculates the potential induced assuming a continuous distribution of doublet strength.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        phi_d : float
            Doublet-induced potential.
        """

        # Calculate geometry
        P = np.array(P)
        geom = self.calc_geom(P)

        # Calculate necessary ints
        I = self.calc_F_integrals(geom)
        self.calc_H_integrals(geom, I)

        # Calculate potential
        phi_d = (self.mu_params[0]*I.hH113 # mu_0
                 + self.mu_params[1]*(P[0]*I.hH113 + geom.h*I.H213) # mu_x
                 + self.mu_params[2]*(P[1]*I.hH113 + geom.h*I.H123) # mu_y
                 + self.mu_params[3]*(0.5*P[0]**2*I.hH113 + geom.h*(P[0]*I.H213 + 0.5*I.H313)) # mu_xx
                 + self.mu_params[4]*(P[0]*P[1]*I.hH113 + geom.h*(P[0]*I.H123 + P[1]*I.H213 + I.H223)) # mu_xy
                 + self.mu_params[5]*(0.5*P[0]**2*I.hH113 + geom.h*(P[1]*I.H123 + 0.5*I.H133)) # mu_yy
                 )*K_inv

        return phi_d


    def calc_discrete_source_potential(self, P, Nx, Ny):
        """Calculates the potential induced assuming a distribution of discrete sources across the panel surface.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.
            
        Nx : integer
            Number of point sources to distribute in the x direction.
            
        Ny : integer
            Number of point sources to distribute in the y direction.

        Returns
        -------
        phi_s : float
            source-induced potential.
        """

        # Distribute sources
        N = Nx*Ny
        X = np.linspace(-0.5*self.x_dim, 0.5*self.x_dim, Nx)
        Y = np.linspace(-0.5*self.y_dim, 0.5*self.y_dim, Ny)

        # Loop through sources
        phi_s = 0.0
        for i, xi in enumerate(X):
            for j, yj in enumerate(Y):

                # Get source strength
                sigma = self.sigma_params[0] + self.sigma_params[1]*xi + self.sigma_params[2]*yj

                # Calculate induced potential
                R = np.sqrt((P[0]-xi)**2 + (P[1]-yj)**2 + P[2]**2)
                phi_s += -sigma*K_inv/R

        return phi_s/N


    def calc_discrete_doublet_potential(self, P, Nx, Ny):
        """Calculates the potential induced assuming a distribution of discrete doublets across the panel surface.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.
            
        Nx : integer
            Number of point doublets to distribute in the x direction.
            
        Ny : integer
            Number of point doublets to distribute in the y direction.

        Returns
        -------
        phi_d : float
            Doublet-induced potential.
        """

        # Distribute doublets
        N = Nx*Ny
        X = np.linspace(-0.5*self.x_dim, 0.5*self.x_dim, Nx)
        Y = np.linspace(-0.5*self.y_dim, 0.5*self.y_dim, Ny)

        # Loop through doublets
        phi_d = 0.0
        for i, xi in enumerate(X):
            for j, yj in enumerate(Y):

                # Get doublet strength
                mu = self.mu_params[0] + self.mu_params[1]*xi + self.mu_params[2]*yj + 0.5*self.mu_params[3]*xi**2 + self.mu_params[4]*xi*yj + 0.5*self.mu_params[5]*yj**2

                # Calculate induced potential
                R = np.sqrt((P[0]-xi)**2 + (P[1]-yj)**2 + P[2]**2)
                phi_d += mu*P[2]*K_inv/R**3

        return phi_d/N


class SupersonicSubinclinedPanel:
    """A class defining a rectangular, supersonic, subinclined, quadratic-doublet-linear-source panel in incompressible flow.
    The panel lies in the x-y (z=0) plane and is centered at the origin. The freestream is aligned with the x-axis.

           0
       ----------
       |        |
     1 |        | 3
       |        |
       ----------
            2
    
    Parameters
    ----------
    x_dim : float
        Width in the x direction.

    y_dim : float
        Width in the y direction.
    """

    def __init__(self, x_dim, y_dim):

        # Store
        self.x_dim = x_dim
        self.y_dim = y_dim
        self.verts = np.zeros((2,4))
        
        # First vertex
        self.verts[0,0] = 0.5*x_dim
        self.verts[1,0] = 0.5*y_dim
        
        # Second vertex
        self.verts[0,1] = -0.5*x_dim
        self.verts[1,1] = 0.5*y_dim
        
        # Third vertex
        self.verts[0,2] = -0.5*x_dim
        self.verts[1,2] = -0.5*y_dim
        
        # Fourth vertex
        self.verts[0,3] = 0.5*x_dim
        self.verts[1,3] = -0.5*y_dim

        # Initialize a few things
        self.mu_params = np.zeros(6)
        self.sigma_params = np.zeros(3)


    def set_doublet_strength(self, mu_params):
        self.mu_params = mu_params


    def set_source_strength(self, sigma_params):
        self.sigma_params = sigma_params


    def calc_geom(self, P):
        """Calculates the geometry for the panel-point pair.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        geom : Geometry
            Geometry of the point relative to the panel
        """

        # Initialize
        geom = Geometry()
        geom.h = P[2]
        geom.h2 = geom.h**2

        # Outward normals
        geom.v_eta[0] = 1.0
        geom.v_eta[2] = -1.0
        geom.v_xi[1] = 1.0
        geom.v_xi[3] = -1.0

        # Edge inclination
        geom.b = geom.v_xi**2 - geom.v_eta**2

        # Loop through edges
        for i in range(4):

            # Displacement to first vertex
            d = self.verts[:,i] - P[0:2]
            
            # Perpendicular in-plane distance
            geom.a[i] = d[0]*geom.v_xi[i] + d[1]*geom.v_eta[i]

            # Tangential in-plane distance
            geom.l1[i] = d[0]*geom.v_eta[i] + d[1]*geom.v_xi[i]

            # Perpendicular distance
            geom.g2[i] = geom.a[i]**2 - geom.b[i]*geom.h2

            # Radial distance
            R1_sq = d[0]**2 - d[1]**2 - geom.h2
            if R1_sq > 0.0 and d[0] < 0.0:
                geom.R1[i] = np.sqrt(R1_sq)
            else:
                geom.R1[i] = 0.0
                geom.l1[i] = -np.sqrt(np.abs(geom.g2[i]))

            # Displacement to second vertex
            i_next = (i+1)%4
            d = self.verts[:,i_next] - P[0:2]

            # Tangential in-plane distance
            geom.l2[i] = d[0]*geom.v_eta[i] + d[1]*geom.v_xi[i]

            # Radial distance
            R2_sq = d[0]**2 - d[1]**2 - geom.h2
            if R2_sq > 0.0 and d[0] < 0.0:
                geom.R2[i] = np.sqrt(R2_sq)
            else:
                geom.R2[i] = 0.0
                geom.l2[i] = np.sqrt(np.abs(geom.g2[i]))

        return geom


    def calc_F_integrals(self, geom):
        """Calculates necessary F integrals.
        
        Parameters
        ----------
        geom : Geometry
            Geometry of the point relative to the panel

        Returns
        -------
        ints : Integrals
            Container of F ints.
        """

        ints = Integrals()

        # Loop through edges
        for i in range(4):

            # Check DoD
            if geom.R1[i] > 0.0 and geom.R2[i] > 0.0:

                # Get square root of b
                sb = np.sqrt(np.abs(geom.b[i]))

                # Calculate F factors
                if geom.b[i] > 0.0:
                    F1 = (geom.l1[i]*geom.R2[i] - geom.l2[i]*geom.R1[i]) / geom.g2[i]
                    F2 = (geom.b[i]*geom.R1[i]*geom.R2[i] + geom.l1[i]*geom.l2[i]) / geom.g2[i]
                else:
                    F1 = (geom.R2[i] - geom.R1[i])*(geom.R2[i] + geom.R1[i]) / (geom.l1[i]*geom.R2[i] + geom.l2[i]*geom.R1[i])
                    F2 = (geom.g2[i] - geom.l1[i]**2 - geom.l2[i]**2) / (geom.b[i]*geom.R1[i]*geom.R2[i] - geom.l1[i]*geom.l2[i])

                # Calculate F(1,1,1)

                # Nearly-sonic edge
                if np.abs(F2) > 100.0*np.abs(sb*F1):

                    # Series solution
                    eps = F1/F2
                    eps2 = eps*eps
                    series = eps*eps2*(1.0/3.0 - geom.b*eps2/5.0 + (geom.b[i]*eps2)*(geom.b*eps2)/7.0)
                    ints.F111[i] = -eps + geom.b[i]*series

                # Supersonic edge
                elif geom.b[i] > 0.0:

                    # Mach wedge region
                    if geom.R1[i] == 0.0 and geom.R2[i] == 0.0:
                        ints.F111[i] = np.pi/sb

                    # At least one endpoint in
                    else:
                        ints.F111[i] = -np.arctan2(sb*F1, F2)/sb

                # Subsonic edge
                else:
                    F1 = sb*geom.R1[i] + np.abs(geom.l1[i])
                    F2 = sb*geom.R2[i] + np.abs(geom.l2[i])
                    ints.F111[i] = -np.sign(geom.v_eta[i])*np.log(F1/F2)

        # Calculate F(2,1,1) and F(1,2,1)
        ints.F121 = -(geom.a*geom.v_eta*ints.F111 + geom.v_xi*(geom.R2-geom.R1))/geom.b
        ints.F211 = geom.a*geom.v_xi*ints.F111 - geom.v_eta*(geom.R2-geom.R1)

        return ints


    def calc_H_integrals(self, geom, ints):
        """Calculates H(M,N,K) integrals.
        
        Parameters
        ----------
        geom : Geometry
            Geometry of the point relative to the panel

        ints : Integrals
            Integrals calculated thus far
        """

        # Loop through edges
        ints.hH113 = 0.0
        for i in range(4):

            # Check DoD
            if geom.R1[i] > 0.0 and geom.R2[i] > 0.0:

                # Get square root of b
                sb = np.sqrt(np.abs(geom.b[i]))

            # Check in plane of panel
            if abs(geom.h) > 1e-12:

                # F factors
                if geom.b[i] > 0.0:
                    F1 = (geom.l1[i]*geom.R2[i] - geom.l2[i]*geom.R1[i]) / geom.g2[i]
                    F2 = (geom.b[i]*geom.R1[i]*geom.R2[i] + geom.l1[i]*geom.l2[i]) / geom.g2[i]
                else:
                    F1 = (geom.R2[i] - geom.R1[i])*(geom.R2[i] + geom.R1[i]) / (geom.l1[i]*geom.R2[i] + geom.l2[i]*geom.R1[i])
                    F2 = (geom.g2[i] - geom.l1[i]**2 - geom.l2[i]**2) / (geom.b[i]*geom.R1[i]*geom.R2[i] - geom.l1[i]*geom.l2[i])

                # Supersonic edge
                if geom.b[i] > 0.0:

                    # Mach wedge region
                    if geom.R1[i] == 0.0 and geom.R2[i] == 0.0:
                        ints.hH113 += np.pi*np.sign(geom.h*geom.v_xi[i])

                    # Otherwise
                    else:
                        ints.hH113 += np.arctan2(geom.h*geom.a[i]*F1, geom.R1[i]*geom.R2[i] + geom.h2*F2)

                # Subsonic edge
                else:
                    ints.hH113 += np.arctan2(geom.h*geom.a[i]*F1, geom.R1[i]*geom.R2[i] + geom.h2*F2)

        # Calculate H(1,1,1)
        ints.H111 = geom.h*ints.hH113 + np.sum(geom.a*ints.F111).item()

        # Calcualte H(2,1,3) and H(1,2,3)
        ints.H213 = -np.sum(geom.v_xi*ints.F111).item()
        ints.H123 = np.sum(geom.v_eta*ints.F111).item()
        
        # Calculate H(2,1,1) and H(1,2,1)
        ints.H211 = 0.5*(geom.h2*ints.H213 + np.sum(geom.a*ints.F211).item())
        ints.H121 = 0.5*(geom.h2*ints.H123 + np.sum(geom.a*ints.F121).item())

        # Calculate H(3,1,3), H(2,2,3), and H(1,3,3)
        ints.H313 = ints.H111 - np.sum(geom.v_xi*ints.F211).item()
        ints.H223 = -np.sum(geom.v_xi*ints.F121).item()
        ints.H133 = -ints.H111 + np.sum(geom.v_eta*ints.F121).item()

        # Check based on (E5) and (E6)
        assert(abs(-np.sum(geom.v_xi*ints.F121).item() - np.sum(geom.v_eta*ints.F211).item()) < 1e-12)

        # Check based on (E4)
        assert(abs(ints.H111 - ints.H313 + ints.H133 + geom.h*ints.hH113) < 1e-12)


    def calc_analytic_source_potential(self, P):
        """Calculates the potential induced assuming a continuous distribution of source strength.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        phi_s : float
            source-induced potential.
        """

        # Calculate geometry
        P = np.array(P)
        geom = self.calc_geom(P)

        # Calculate necessary ints
        I = self.calc_F_integrals(geom)
        self.calc_H_integrals(geom, I)

        phi_s = self.sigma_params[0]*I.H111 + self.sigma_params[1]*(I.H111*P[0] + I.H211) + self.sigma_params[2]*(I.H111*P[1] + I.H121)
        return -phi_s*K_inv*2.0


    def calc_analytic_doublet_potential(self, P):
        """Calculates the potential induced assuming a continuous distribution of doublet strength.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        phi_d : float
            Doublet-induced potential.
        """

        # Calculate geometry
        P = np.array(P)
        geom = self.calc_geom(P)

        # Calculate necessary ints
        I = self.calc_F_integrals(geom)
        self.calc_H_integrals(geom, I)

        # Calculate potential
        phi_d = ( self.mu_params[0]*I.hH113 # mu_0
                + self.mu_params[1]*(P[0]*I.hH113 + geom.h*I.H213) # mu_x
                + self.mu_params[2]*(P[1]*I.hH113 + geom.h*I.H123) # mu_y
                + self.mu_params[3]*(0.5*P[0]**2*I.hH113 + geom.h*(P[0]*I.H213 + 0.5*I.H313)) # mu_xx
                + self.mu_params[4]*(P[0]*P[1]*I.hH113 + geom.h*(P[0]*I.H123 + P[1]*I.H213 + I.H223)) # mu_xy
                + self.mu_params[5]*(0.5*P[1]**2*I.hH113 + geom.h*(P[1]*I.H123 + 0.5*I.H133)) ) # mu_yy

        return phi_d*K_inv*2.0


    def calc_discrete_source_potential(self, P, Nx, Ny):
        """Calculates the potential induced assuming a distribution of discrete sources across the panel surface.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.
            
        Nx : integer
            Number of point sources to distribute in the x direction.
            
        Ny : integer
            Number of point sources to distribute in the y direction.

        Returns
        -------
        phi_s : float
            source-induced potential.
        """

        # Distribute sources
        N = Nx*Ny
        X = np.linspace(-0.5*self.x_dim, 0.5*self.x_dim, Nx)
        Y = np.linspace(-0.5*self.y_dim, 0.5*self.y_dim, Ny)

        # Loop through sources
        phi_s = 0.0
        for i, xi in enumerate(X):
            for j, yj in enumerate(Y):

                # Calculate hyperbolic distance
                R2 = (P[0]-xi)**2 - (P[1]-yj)**2 - P[2]**2

                # Check
                if R2 > 0.0 and P[0]-xi > 0.0:

                    # Get distance
                    R = np.sqrt(R2)

                    # Get source strength
                    sigma = self.sigma_params[0] + self.sigma_params[1]*xi + self.sigma_params[2]*yj

                    # Calculate induced potential
                    phi_s += sigma/R

        return -2.0*K_inv*phi_s/N


    def calc_discrete_doublet_potential(self, P, Nx, Ny):
        """Calculates the potential induced assuming a distribution of discrete doublets across the panel surface.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.
            
        Nx : integer
            Number of point doublets to distribute in the x direction.
            
        Ny : integer
            Number of point doublets to distribute in the y direction.

        Returns
        -------
        phi_d : float
            Doublet-induced potential.
        """

        # Distribute doublets
        N = Nx*Ny
        X = np.linspace(-0.5*self.x_dim, 0.5*self.x_dim, Nx)
        Y = np.linspace(-0.5*self.y_dim, 0.5*self.y_dim, Ny)

        # Loop through doublets
        phi_d = 0.0
        for i, xi in enumerate(X):
            for j, yj in enumerate(Y):

                # Calculate hyperbolic distance
                R2 = (P[0]-xi)**2 - (P[1]-yj)**2 - P[2]**2

                # Check
                if R2 > 0.0:

                    # Get distance
                    R = np.sqrt(R2)

                    # Get doublet strength
                    mu = self.mu_params[0] + self.mu_params[1]*xi + self.mu_params[2]*yj + 0.5*self.mu_params[3]*xi**2 + self.mu_params[4]*xi*yj + 0.5*self.mu_params[5]*yj**2

                    # Calculate induced potential
                    phi_d += mu/R**3

        return 2.0*K_inv*P[2]*phi_d/N