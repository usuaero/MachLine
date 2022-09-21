import numpy as np
import matplotlib.pyplot as plt


def inner2(x,y):
    return x[0]*y[0] + x[1]*y[1]


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
        self.t = np.zeros((2,4))
        self.h = 0.0
        self.h2 = 0.0
        self.b = np.ones(4)
        self.edge_in = np.zeros(4, dtype=bool)


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


class Panel:
    """A class defining a quadrlateral panel.
    
    Parameters
    ----------
    verts : ndarray
        Array of vertex locations. Of shape (2,4).
    """

    
    def __init__(self, verts):
        
        # Store
        self.verts = verts

        # Calculate area
        self.A = 0.5*( (self.verts[0,1]-self.verts[0,0])*(self.verts[1,2]-self.verts[1,1])
                      -(self.verts[1,1]-self.verts[1,0])*(self.verts[0,2]-self.verts[0,1])
                      +(self.verts[0,3]-self.verts[0,2])*(self.verts[1,0]-self.verts[1,3])
                      -(self.verts[1,3]-self.verts[1,2])*(self.verts[0,0]-self.verts[0,3]) )


    def get_point_on_surface(self, u, v):
        """Returns the coordinates of a point having the nondimensional coordinates (u, v).
        u and v vary from 0 to 1.

        Parameters
        ----------
        u, v : float
            Nondimensional coordinates.

        Returns
        -------
        point : ndarray
            Surface coordinates of the point
        """

        # Calculate weights
        w = np.zeros(4)
        w[0] = (1.0-u)*(1.0-v)
        w[1] = (1.0-u)*v
        w[2] = u*v
        w[3] = u*(1.0-v)

        # Calcluate point
        point = np.zeros(2)
        point[0] = np.inner(self.verts[0,:], w)
        point[1] = np.inner(self.verts[1,:], w)

        return point


    def get_points_dist(self, Nu, Nv):
        """Returns a set of points distributed across the panel.
        
        Parameters
        ----------
        Nu : int
            Number of points in the first dimension.
            
        Nv : int
            Number of points in the second dimension.
        
        Returns
        -------
        points :: ndarray
            Array of points with shape (N1*N2,2).
        """

        # Get coordinate distributions
        us = np.linspace(0.0, 1.0, Nu)
        vs = np.linspace(0.0, 1.0, Nv)
        Us, Vs = np.meshgrid(us, vs)
        points = np.zeros((Nu*Nv,2))

        # Get points
        for i, (u, v) in enumerate(zip(Us.flatten(), Vs.flatten())):
            points[i,:] = self.get_point_on_surface(u, v)

        return points


    def distribute_points(self, Nu, Nv):
        """Creates a distribution of discrete points across the panel surface.
        
        Parameters
        ----------
        Nu : int
            Number of points in the first dimension.
            
        Nv : int
            Number of points in the second dimension.
        """

        self.points = self.get_points_dist(Nu, Nv)
        self.N = Nu*Nv
        self.dA = self.A/self.N


    def get_local_doublet_strength(self, point):
        """Returns the local doublet strength at the specified point.
        
        Parameters
        ----------
        point : array
            Coordinates of point on surface.

        Returns
        -------
        mu : float
            Doublet strength
        """

        return self.mu_params[0] + self.mu_params[1]*point[0] + self.mu_params[2]*point[1] + 0.5*self.mu_params[3]*point[0]**2 + self.mu_params[4]*point[0]*point[1] + 0.5*self.mu_params[5]*point[1]**2


    def get_local_source_strength(self, point):
        """Returns the local source strength at the specified point.
        
        Parameters
        ----------
        point : array
            Coordinates of point on surface.

        Returns
        -------
        mu : float
            Source strength
        """

        return self.sigma_params[0] + self.sigma_params[1]*point[0] + self.sigma_params[2]*point[1]


class SubsonicPanel(Panel):
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
        verts = np.zeros((2,4))
        
        # First vertex
        verts[0,0] = 0.5*x_dim
        verts[1,0] = 0.5*y_dim
        
        # Second vertex
        verts[0,1] = -0.5*x_dim
        verts[1,1] = 0.5*y_dim
        
        # Third vertex
        verts[0,2] = -0.5*x_dim
        verts[1,2] = -0.5*y_dim
        
        # Fourth vertex
        verts[0,3] = 0.5*x_dim
        verts[1,3] = -0.5*y_dim
        super().__init__(verts)

        # Initialize a few things
        self.mu_params = np.zeros(6)
        self.sigma_params = np.zeros(3)
        self.K_inv = 0.25/np.pi


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
        geom.P = P

        # Loop through edges
        for i in range(4):
            
            # Calculate outward normals
            t = self.verts[:,(i+1)%4] - self.verts[:,i]
            t = t/np.linalg.norm(t)
            geom.v_xi[i] = t[1]
            geom.v_eta[i] = -t[0]

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
        return -phi_s*self.K_inv


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
                 )*self.K_inv

        return phi_d


    def calc_discrete_source_potential(self, P):
        """Calculates the potential induced assuming a distribution of discrete sources across the panel surface.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        phi_s : float
            source-induced potential.
        """

        # Loop through sources
        phi_s = 0.0
        for point in self.points:

            # Calculate induced potential
            R = np.sqrt((P[0]-point[0])**2 + (P[1]-point[1])**2 + P[2]**2)
            phi_s += -self.get_local_source_strength(point)/R

        return phi_s*self.dA*self.K_inv


    def calc_discrete_doublet_potential(self, P):
        """Calculates the potential induced assuming a distribution of discrete doublets across the panel surface.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        phi_d : float
            Doublet-induced potential.
        """

        # Loop through doublets
        phi_d = 0.0
        for point in self.points:

            # Calculate induced potential
            R = np.sqrt((P[0]-point[0])**2 + (P[1]-point[1])**2 + P[2]**2)
            phi_d += self.get_local_doublet_strength(point)/R**3

        return phi_d*P[2]*self.K_inv*self.dA


class SupersonicSubinclinedPanel(Panel):
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
    verts : ndarray
        Array of vertex locations. Of shape (2,4).
    """

    def __init__(self, verts):

        # Store
        super().__init__(verts)

        # Initialize a few things
        self.mu_params = np.zeros(6)
        self.sigma_params = np.zeros(3)
        self.K_inv = 0.5/np.pi


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
        geom.P = P

        # Loop through edges
        for i in range(4):
            
            # Calculate outward normals
            t = self.verts[:,(i+1)%4] - self.verts[:,i]
            t = t/np.linalg.norm(t)
            geom.v_xi[i] = t[1]
            geom.v_eta[i] = -t[0]

            # Inclination parameter
            geom.b[i] = geom.v_xi[i]**2 - geom.v_eta[i]**2

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
            if (geom.R1[i] > 0.0 or geom.R2[i] > 0.0) or (geom.b[i] > 0.0 and geom.g2[i] > 0.0 and geom.l1[i]*geom.l2[i] < 0.0 and geom.a[i]*geom.v_xi[i] < 0.0):

                # Get square root of b
                sb = np.sqrt(np.abs(geom.b[i]))

                # Mach wedge region
                if geom.R1[i] == 0.0 and geom.R2[i] == 0.0:
                    ints.F111[i] = np.pi/sb

                    # F(1,2,1) and F(2,1,1)
                    ints.F121[i] = -geom.a[i]*geom.v_eta[i]*ints.F111[i]/geom.b[i]
                    ints.F211[i] = geom.a[i]*geom.v_xi[i]*ints.F111[i]/geom.b[i]

                else:

                    # Calculate F factors
                    if geom.b[i] > 0.0:
                        F1 = (geom.l1[i]*geom.R2[i] - geom.l2[i]*geom.R1[i]) / geom.g2[i]
                        F2 = (geom.b[i]*geom.R1[i]*geom.R2[i] + geom.l1[i]*geom.l2[i]) / geom.g2[i]
                    else:
                        F1 = (geom.R2[i] - geom.R1[i])*(geom.R2[i] + geom.R1[i]) / (geom.l1[i]*geom.R2[i] + geom.l2[i]*geom.R1[i])
                        F2 = (geom.g2[i] - geom.l1[i]**2 - geom.l2[i]**2) / (geom.b[i]*geom.R1[i]*geom.R2[i] - geom.l1[i]*geom.l2[i])

                    # Calculate F(1,1,1)

                    # Nearly-sonic edge
                    if abs(F2) > 100.0*sb*abs(F1):

                        # Series solution
                        eps = F1/F2
                        eps2 = eps*eps
                        series = eps*eps2*(1.0/3.0 - geom.b[i]*eps2/5.0 + (geom.b[i]*eps2)*(geom.b[i]*eps2)/7.0)
                        ints.F111[i] = -eps + geom.b[i]*series

                        # F(1,2,1) and F(2,1,1)
                        ints.F121[i] = (- geom.v_xi[i]*(geom.R2[i]-geom.R1[i])*geom.R1[i]*geom.R2[i] 
                                        + geom.l2[i]*geom.R1[i]*(self.verts[1,i]-geom.P[1]) 
                                        - geom.l1[i]*geom.R2[i]*(self.verts[1,(i+1)%4]-geom.P[1]) ) / (geom.g2[i]*F2) \
                                        - geom.a[i]*geom.v_eta[i]*series
                        ints.F211 = geom.a*geom.v_xi*ints.F111 - geom.v_eta*(geom.R2-geom.R1) - 2.0*geom.v_xi*geom.v_eta*ints.F121

                    # Supersonic edge
                    elif geom.b[i] > 0.0:
                        ints.F111[i] = -np.arctan2(sb*F1, F2)/sb

                        # F(1,2,1) and F(2,1,1)
                        ints.F121[i] = -(geom.v_xi[i]*(geom.R2[i]-geom.R1[i]) + geom.a[i]*geom.v_eta[i]*ints.F111[i])/geom.b[i]
                        ints.F211 = geom.a*geom.v_xi*ints.F111 - geom.v_eta*(geom.R2-geom.R1) - 2.0*geom.v_xi*geom.v_eta*ints.F121

                    # Subsonic edge
                    else:
                        F1 = sb*geom.R1[i] + np.abs(geom.l1[i])
                        F2 = sb*geom.R2[i] + np.abs(geom.l2[i])
                        ints.F111[i] = -np.sign(geom.v_eta[i])*np.log(F1/F2)/sb

                        # F(1,2,1) and F(2,1,1)
                        ints.F121[i] = -(geom.v_xi[i]*(geom.R2[i]-geom.R1[i]) + geom.a[i]*geom.v_eta[i]*ints.F111[i])/geom.b[i]
                        ints.F211 = geom.a*geom.v_xi*ints.F111 - geom.v_eta*(geom.R2-geom.R1) - 2.0*geom.v_xi*geom.v_eta*ints.F121
        
        # Check
        assert((np.abs(geom.v_xi*ints.F211 + geom.v_eta*ints.F121 - geom.a*ints.F111) < 1.0e-12).all())

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
            if (geom.R1[i] > 0.0 or geom.R2[i] > 0.0) or (geom.b[i] > 0.0 and geom.g2[i] > 0.0 and geom.l1[i]*geom.l2[i] < 0.0 and geom.a[i]*geom.v_xi[i] < 0.0):

                # Check in plane of panel
                if abs(geom.h) > 1e-12:

                    # Mach wedge region
                    if geom.R1[i] == 0.0 and geom.R2[i] == 0.0:
                        ints.hH113 -= np.pi*np.sign(geom.h*geom.v_xi[i])

                    else:

                        # F factors
                        if geom.b[i] > 0.0:
                            F1 = (geom.l1[i]*geom.R2[i] - geom.l2[i]*geom.R1[i]) / geom.g2[i]
                            F2 = (geom.b[i]*geom.R1[i]*geom.R2[i] + geom.l1[i]*geom.l2[i]) / geom.g2[i]
                        else:
                            F1 = (geom.R2[i] - geom.R1[i])*(geom.R2[i] + geom.R1[i]) / (geom.l1[i]*geom.R2[i] + geom.l2[i]*geom.R1[i])
                            F2 = (geom.g2[i] - geom.l1[i]**2 - geom.l2[i]**2) / (geom.b[i]*geom.R1[i]*geom.R2[i] - geom.l1[i]*geom.l2[i])

                        ints.hH113 -= np.arctan2(geom.h*geom.a[i]*F1, geom.R1[i]*geom.R2[i] + geom.h2*F2)

        # Calculate H(1,1,1)
        ints.H111 = -geom.h*ints.hH113 + np.sum(geom.a*ints.F111).item()

        # Calcualte H(2,1,3) and H(1,2,3)
        ints.H213 = np.sum(geom.v_xi*ints.F111).item()
        ints.H123 = -np.sum(geom.v_eta*ints.F111).item()
        
        # Calculate H(2,1,1) and H(1,2,1)
        ints.H211 = 0.5*(-geom.h2*ints.H213 + np.sum(geom.a*ints.F211).item())
        ints.H121 = 0.5*(-geom.h2*ints.H123 + np.sum(geom.a*ints.F121).item())

        # Calculate H(3,1,3), H(2,2,3), and H(1,3,3)
        ints.H313 = -ints.H111 + np.sum(geom.v_xi*ints.F211).item()
        ints.H223 = np.sum(geom.v_xi*ints.F121).item()
        ints.H133 = ints.H111 - np.sum(geom.v_eta*ints.F121).item()

        # Check based on (E5) and (E6)
        assert(abs(ints.H223 + np.sum(geom.v_eta*ints.F211).item()) < 1e-12)

        # Check based on (E4)
        assert(abs(ints.H111 + ints.H313 - ints.H133 - geom.h*ints.hH113) < 1e-12)


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
        return -phi_s*self.K_inv*2.0


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

        return phi_d*self.K_inv*2.0


    def calc_discrete_source_potential(self, P):
        """Calculates the potential induced assuming a distribution of discrete sources across the panel surface.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        phi_s : float
            source-induced potential.
        """

        # Loop through sources
        phi_s = 0.0
        for point in self.points:

            # Calculate hyperbolic distance
            R2 = (P[0]-point[0])**2 - (P[1]-point[1])**2 - P[2]**2

            # Check
            if R2 > 0.0 and P[0]-point[0] > 0.0:

                # Get distance
                R = np.sqrt(R2)

                # Calculate induced potential
                phi_s += self.get_local_source_strength(point)/R

        return -2.0*self.K_inv*phi_s*self.dA


    def calc_discrete_doublet_potential(self, P):
        """Calculates the potential induced assuming a distribution of discrete doublets across the panel surface.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        phi_d : float
            Doublet-induced potential.
        """

        # Loop through doublets
        phi_d = 0.0
        for point in self.points:

            # Calculate hyperbolic distance
            R2 = (P[0]-point[0])**2 - (P[1]-point[1])**2 - P[2]**2

            # Check
            if R2 > 0.0 and P[0]-point[0] > 0.0:

                # Get distance
                R = np.sqrt(R2)

                # Calculate induced potential
                phi_d += self.get_local_doublet_strength(point)/R**3

        return -2.0*self.K_inv*P[2]*phi_d*self.dA # I don't know why a negative sign is needed here, but it makes everything in MachLine work


class SupersonicSuperinclinedPanel(Panel):
    """A class defining a rectangular, supersonic, subinclined, quadratic-doublet-linear-source panel in incompressible flow.
    The panel lies in the x-y (z=0) plane and is centered at the origin. The freestream is aligned with the z-axis.

           0
       ----------
       |        |
     1 |        | 3
       |        |
       ----------
            2
    
    Parameters
    ----------
    verts : ndarray
        Array of vertex locations. Of shape (2,4).
    """

    def __init__(self, verts):

        # Store
        super().__init__(verts)

        # Initialize a few things
        self.mu_params = np.zeros(6)
        self.sigma_params = np.zeros(3)
        self.K_inv = 0.5/np.pi


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
        geom.P = P
        geom.b = 1.0 # Meaningless really, but we'll keep it here

        # Loop through edges
        for i in range(4):
            
            # Calculate tangents
            i_next = (i+1)%4
            geom.t[:,i] = self.verts[:,i_next] - self.verts[:,i]
            geom.t[:,i] = geom.t[:,i]/np.linalg.norm(geom.t[:,i])

            # Calculate outward normals
            geom.v_xi[i] = geom.t[1,i]
            geom.v_eta[i] = -geom.t[0,i]

            # Displacement to first vertex
            d = self.verts[:,i] - P[0:2]
            
            # Perpendicular in-plane distance
            geom.a[i] = d[0]*geom.v_xi[i] + d[1]*geom.v_eta[i]

            # Tangential in-plane distance
            geom.l1[i] = d[0]*geom.t[0,i] + d[1]*geom.t[1,i]

            # Radial distance
            R1_sq = -d[0]**2 - d[1]**2 + geom.h2

            # Displacement to second vertex
            d = self.verts[:,i_next] - P[0:2]

            # Tangential in-plane distance
            geom.l2[i] = d[0]*geom.t[0,i] + d[1]*geom.t[1,i]

            # Radial distance
            R2_sq = -d[0]**2 - d[1]**2 + geom.h2

            # Check if edge is in
            if geom.h > 0.0 and (R1_sq > 0.0 or R2_sq > 0.0) or (abs(geom.a[i]) < geom.h and geom.l1[i]*geom.l2[i] < 0.0):
                geom.edge_in[i] = True
            else:
                geom.edge_in[i] = False

            # Check corners for DoD
            if R1_sq > 0.0 and geom.h > 0.0:
                geom.R1[i] = np.sqrt(R1_sq)
            else:
                geom.R1[i] = 0.0
                geom.l1[i] = -1.0

            # Check corners for DoD
            if R2_sq > 0.0 and geom.h > 0.0:
                geom.R2[i] = np.sqrt(R2_sq)
            else:
                geom.R2[i] = 0.0
                geom.l2[i] = 1.0

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
            if geom.edge_in[i]:

                # Mach wedge
                if geom.R1[i] == 0.0 and geom.R2[i] == 0.0:
                    ints.F111[i] = -np.pi

                    # F(1,2,1) and F(2,1,1)
                    ints.F211[i] = geom.a[i]*geom.v_xi[i]*ints.F111[i]
                    ints.F121[i] = geom.a[i]*geom.v_eta[i]*ints.F111[i]

                # At least one endpoint in
                else:
                    X = geom.l1[i]*geom.l2[i] + geom.R1[i]*geom.R2[i]
                    Y = geom.R2[i]*geom.l1[i] - geom.R1[i]*geom.l2[i]
                    ints.F111[i] = np.arctan2(Y, X)

                    # F(1,2,1) and F(2,1,1)
                    ints.F211[i] = geom.a[i]*geom.v_xi[i]*ints.F111[i] - geom.v_eta[i]*(geom.R2[i]-geom.R1[i])
                    ints.F121[i] = geom.a[i]*geom.v_eta[i]*ints.F111[i] + geom.v_xi[i]*(geom.R2[i]-geom.R1[i])
        # Check
        assert((np.abs(geom.v_xi*ints.F211 + geom.v_eta*ints.F121 - geom.a*ints.F111) < 1.0e-12).all())

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

        # Check DoD
        if geom.edge_in.any() or (not geom.edge_in.any() and (geom.a > 0.0).all() and geom.h > 0.0):

            # Loop through edges
            ints.hH113 = 2.0*np.pi
            for i in range(4):

                # Add corner influence
                if geom.R1[i] > 0.0:
                    ints.hH113 += np.pi

                # Add edge influence
                if geom.edge_in[i]:
                    ints.hH113 -= np.pi

                # Add complicated corner influence
                if geom.R1[i] > 0.0:
                    i_prev = i-1

                    # Get dot and cross products
                    tktkp1 = np.inner(geom.t[:,i], geom.t[:,i_prev])
                    tcross = np.cross(geom.t[:,i_prev], geom.t[:,i])

                    # Intermediate values
                    X = geom.a[i]*geom.a[i_prev] - geom.h2*tktkp1
                    Y = geom.h*geom.R1[i]*tcross

                    # Update hH113
                    ints.hH113 -= np.arctan2(Y, -X)

        # Calculate H(1,1,1)
        ints.H111 = geom.h*ints.hH113 - np.sum(geom.a*ints.F111).item()

        # Calcualte H(2,1,3) and H(1,2,3)
        ints.H213 = np.sum(geom.v_xi*ints.F111).item()
        ints.H123 = np.sum(geom.v_eta*ints.F111).item()
        
        # Calculate H(2,1,1) and H(1,2,1)
        ints.H211 = 0.5*(geom.h2*ints.H213 - np.sum(geom.a*ints.F211).item())
        ints.H121 = 0.5*(geom.h2*ints.H123 - np.sum(geom.a*ints.F121).item())

        # Calculate H(3,1,3), H(2,2,3), and H(1,3,3)
        ints.H313 = ints.H111 + np.sum(geom.v_xi*ints.F211).item()
        ints.H223 = np.sum(geom.v_xi*ints.F121).item()
        ints.H133 = ints.H111 + np.sum(geom.v_eta*ints.F121).item()

        # Check based on (E5) and (E6)
        assert(abs(ints.H223 - np.sum(geom.v_eta*ints.F211).item()) < 1e-12)

        # Check based on (E4)
        assert(abs(ints.H111 - ints.H313 - ints.H133 + geom.h*ints.hH113) < 1e-12)


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
        return -phi_s*self.K_inv*2.0


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

        return phi_d*self.K_inv*2.0


    def calc_discrete_source_potential(self, P):
        """Calculates the potential induced assuming a distribution of discrete sources across the panel surface.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        phi_s : float
            source-induced potential.
        """

        # Loop through sources
        phi_s = 0.0
        for point in self.points:

            # Calculate hyperbolic distance
            R2 = -(P[0]-point[0])**2 - (P[1]-point[1])**2 + P[2]**2

            # Check
            if R2 > 0.0 and P[2] > 0.0:

                # Get distance
                R = np.sqrt(R2)

                # Calculate induced potential
                phi_s += self.get_local_source_strength(point)/R

        return -2.0*self.K_inv*phi_s*self.dA


    def calc_discrete_doublet_potential(self, P):
        """Calculates the potential induced assuming a distribution of discrete doublets across the panel surface.
        
        Parameters
        ----------
        P : ndarray
            Point at which to calculate the induced potential.

        Returns
        -------
        phi_d : float
            Doublet-induced potential.
        """

        # Loop through doublets
        phi_d = 0.0
        for point in self.points:

            # Calculate hyperbolic distance
            R2 = -(P[0]-point[0])**2 - (P[1]-point[1])**2 + P[2]**2

            # Check
            if R2 > 0.0 and P[2] > 0.0:

                # Get distance
                R = np.sqrt(R2)

                # Calculate induced potential
                phi_d += self.get_local_doublet_strength(point)/R**3

        return -2.0*self.K_inv*P[2]*phi_d*self.dA # I don't know why a negative sign is needed here, but it makes everything in MachLine work