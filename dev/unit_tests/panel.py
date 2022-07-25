import numpy as np


def inner2(x,y):
    return x[0]*y[0] + x[1]*y[1]


class SubsonicGeometry:
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


class Integrals:
    """A class containing all necessary integrals."""

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
        geom : SubsonicGeometry
            Geometry of the point relative to the panel
        """

        # Initialize
        geom = SubsonicGeometry()
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
        geom : SubsonicGeometry
            Geometry of the point relative to the panel

        Returns
        -------
        integrals : Integrals
            Container of F integrals.
        """

        integrals = Integrals()

        # Loop through edges
        for i in range(4):

            # Within edge
            if geom.l1[i]*geom.l2[i] < 0.0:

                integrals.F111[i] = np.log(((geom.R1[i] - geom.l1[i])*(geom.R2[i] + geom.l2[i]))/geom.g2[i])

            # Outside edge
            else:

                integrals.F111[i] = np.sign(geom.l1[i])*np.log((geom.R2[i]+np.abs(geom.l2[i]))/(geom.R1[i]+np.abs(geom.l1[i])))

        # Calculate F(2,1,1) and F(1,2,1)
        integrals.F121 = geom.a*geom.v_eta*integrals.F111 + geom.v_xi*(geom.R2-geom.R1)
        integrals.F211 = geom.a*geom.v_xi*integrals.F111 - geom.v_eta*(geom.R2-geom.R1)

        return integrals


    def calc_H_integrals(self, geom, integrals):
        """Calculates hH(1,1,3).
        
        Parameters
        ----------
        geom : SubsonicGeometry
            Geometry of the point relative to the panel

        integrals : Integrals
            Integrals calculated thus far
        """

        # Loop through edges
        integrals.hH113 = 0.0
        for i in range(4):
            
            # Intermediate quantities
            c1 = geom.g2[i] + np.abs(geom.h)*geom.R1[i]
            c2 = geom.g2[i] + np.abs(geom.h)*geom.R2[i]

            # Integral for edge
            S = geom.a[i]*(geom.l2[i]*c1 - geom.l1[i]*c2)
            C = c1*c2 + geom.a[i]**2*geom.l1[i]*geom.l2[i]

            # Sum
            integrals.hH113 += np.arctan2(S, C)

        # Apply sign factor
        integrals.hH113 *= np.sign(geom.h)

        # Calculate H(1,1,1)
        integrals.H111 = -geom.h*integrals.hH113 + np.sum(geom.a*integrals.F111).item()

        # Calcualte H(2,1,3) and H(1,2,3)
        integrals.H213 = -np.sum(geom.v_xi*integrals.F111).item()
        integrals.H123 = -np.sum(geom.v_eta*integrals.F111).item()
        
        # Calculate H(2,1,1) and H(1,2,1)
        integrals.H211 = 0.5*(-geom.h**2*integrals.H213 + np.sum(geom.a*integrals.F211).item())
        integrals.H121 = 0.5*(-geom.h**2*integrals.H123 + np.sum(geom.a*integrals.F121).item())


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

        # Calculate necessary integrals
        I = self.calc_F_integrals(geom)
        self.calc_H_integrals(geom, I)

        phi_s = self.sigma_params[0]*I.H111 + self.sigma_params[1]*(I.H111*P[0] + I.H211) + self.sigma_params[2]*(I.H111*P[1] + I.H121)
        return -phi_s/(4.0*np.pi)


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

        # Calculate necessary integrals
        I = self.calc_F_integrals(geom)
        self.calc_H_integrals(geom, I)

        # Calculate potential
        phi_d = (self.mu_params[0]*I.hH113 # mu_0
                 + self.mu_params[1]*(P[0]*I.hH113 + geom.h*I.H213) # mu_x
                 + self.mu_params[2]*(P[1]*I.hH113 + geom.h*I.H123) # mu_y
                 )/(4.0*np.pi)

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
                phi_s += -sigma/(4.0*np.pi*R)

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
                phi_d += mu*P[2]/(4.0*np.pi*R**3)

        return phi_d/N


class SubinclinedPanel:
    """A class defining a subinclined constant-doublet-constant-source panel of unit strength in supersonic flow (M_0=sqrt(2)).
    The flow is aligned with the x axis.
    
    Parameters
    ----------
    verts : ndarray
        2xN array, the columns of which are the corner points of the panel, expressed in local coords.

    mu : float
        Doublet strength.

    sigma : float
        Source strength.
    """

    def __init__(self, verts, mu, sigma):

        # Store
        self.verts = verts
        self.N = self.verts.shape[1]
        self.mu = mu
        self.sigma = sigma

        # Initialize storage
        self.dx = np.zeros(self.N)
        self.dy = np.zeros(self.N)
        self.t_hat = np.zeros((2,self.N))
        self.n_hat = np.zeros((2,self.N))
        self.b = np.zeros(self.N)

        # Calculate edge parameters
        for i in range(self.N):

            # Get index of end vertex
            i_next = (i+1)%self.N

            # Get edge displacements
            self.dx[i] = verts[0,i_next] - verts[0,i]
            self.dy[i] = verts[1,i_next] - verts[1,i]

            # Get edge tangent vector
            self.t_hat[0,i] = self.dx[i]
            self.t_hat[1,i] = self.dy[i]
            self.t_hat[:,i] /= np.linalg.norm(self.t_hat[:,i])

        # Get edge normal vector
        self.n_hat[0,:] = self.t_hat[1,:]
        self.n_hat[1,:] = -self.t_hat[0,:]

        # Get b
        self.b = (self.n_hat[0,:] - self.n_hat[1,:])*(self.n_hat[0,:] + self.n_hat[1,:])

        # Calculate area
        if self.N == 3:
            self.A = 0.5*abs(self.dx[0]*self.dy[1] - self.dy[0]*self.dx[1])
        else:
            self.A = abs(self.dx[0]*self.dy[1] - self.dy[0]*self.dx[1])

        # Calculate centroid
        self.c = np.sum(self.verts, axis=1)/self.N


    def _calc_geometry(self, P):
        # Calculates necessary geometry

        # Set up geometry
        h = P[2]
        R1 = np.zeros(self.N)
        R2 = np.zeros(self.N)
        l1 = np.zeros(self.N)
        l2 = np.zeros(self.N)
        a = np.zeros(self.N)
        g2 = np.zeros(self.N)
        in_dod = np.zeros(self.N, dtype=bool)

        # Loop through edges
        for i in range(self.N):

            # Get index of end vertex
            i_next = (i+1)%self.N

            # Get displacements
            # First vertex
            d1 = self.verts[:,i] - P[0:2]
            l1[i] = self.n_hat[1,i]*d1[0] + self.n_hat[0,i]*d1[1]

            # Second vertex
            d2 = self.verts[:,i_next] - P[0:2]
            l2[i] = self.n_hat[1,i]*d2[0] + self.n_hat[0,i]*d2[1]

            # Perpendicular distances
            a[i] = inner2(self.n_hat[:,i], d1)
            g2[i] = a[i]**2 - self.b[i]*h**2

            # Get hyperbolic radii
            x = d1[0]**2 - d1[1]**2 - h**2
            if x > 0.0 and d1[0] < 0.0:
                R1[i] = np.sqrt(x)
            else:
                R1[i] = 0.0

            x = d2[0]**2 - d2[1]**2 - h**2
            if x > 0.0 and d2[0] < 0.0:
                R2[i] = np.sqrt(x)
            else:
                R2[i] = 0.0

            # Check DoD

            # Check if both endpoints are out
            if R1[i] == 0.0 and R2[i] == 0.0:

                # In this case, a subsonic or sonic edge is out
                if self.b[i] <= 0.0:
                    in_dod[i] = False

                # For a supersonic edge, it may still be in
                else:

                    # Check for outside edge
                    if l1[i]*l2[i] >= 0.0:
                        in_dod[i] = False

                    # Check for above or below Mach wedge
                    elif g2[i] <= 0.0:
                        in_dod[i] = False

                    # Check for upstream
                    elif a[i]*self.n_hat[0,i] >= 0.0:
                        in_dod[i] = False

                    else:
                        in_dod[i] = True

            else:
                in_dod[i] = True

            # Fix l from PAN AIR
            if R1[i] == 0.0:
                l1[i] = -np.sqrt(abs(g2[i]))
            if R2[i] == 0.0:
                l2[i] = np.sqrt(abs(g2[i]))

        return h, R1, R2, l1, l2, a, g2, in_dod


    def _calc_integrals(self, P):
        # Calculates the needed integrals

        # Calculate geometry
        h, R1, R2, l1, l2, a, g2, in_dod = self._calc_geometry(P)

        # Initialize
        hH113 = np.zeros(self.N)
        F111 = np.zeros(self.N)

        # Loop through edges
        for i in range(self.N):

            # Check DoD
            if in_dod[i]:

                # Calculate F factors
                s_b = np.sqrt(abs(self.b[i]))
                if self.b[i] > 0.0:
                    F1 = (l1[i]*R2[i] - l2[i]*R1[i]) / g2[i]
                    F2 = (self.b[i]*R1[i]*R2[i] + l1[i]*l2[i]) / g2[i]
                else:
                    F1 = (R2[i] - R1[i])*(R2[i] + R1[i]) / (l1[i]*R2[i] + l2[i]*R1[i])
                    F2 = (g2[i] - l1[i]**2 - l2[i]**2) / (self.b[i]*R1[i]*R2[i] - l1[i]*l2[i])

                # Calculate hH(1,1,3)

                # Check for point on panel plane
                if h != 0.0:

                    # Supersonic edge
                    if self.b[i] > 0.0:

                        # Neither endpoint in
                        if R1[i] == 0.0 and R2[i] == 0.0:
                            hH113[i] = np.pi*np.sign(h*self.n_hat[0,i])

                        # At least one endpoint in
                        else:
                            hH113[i] = np.arctan2(h*a[i]*F1, R1[i]*R2[i] + h**2*F2)
                
                    # Subsonic edge
                    else:
                        hH113[i] = np.arctan2(h*a[i]*F1, R1[i]*R2[i] + h**2*F2)


                # Calculate F(1,1,1)

                # Nearly-sonic edge
                if abs(F2) > 100.0*s_b*abs(F1):

                    # Use series solution
                    e = F1/F2
                    F111[i] = -e*(1.0 - self.b[i]*e**2/3.0 + self.b[i]**2*e**4/5.0 - self.b[i]**3*e**6/7.0)

                # Supersonic edge
                elif self.b[i] > 0.0:

                    # Check for Mach wedge condition
                    if R1[i] == 0.0 and R2[i] == 0.0:
                        F111[i] = np.pi/s_b

                    # At least one in
                    else:
                        F111[i] = -np.arctan2(s_b*F1, F2) / s_b

                # Subsonic edge
                else:
                    F111[i] = -np.sign(self.n_hat[1,i]) * np.log( (s_b*R1[i] + abs(l1[i])) / (s_b*R2[i] + abs(l2[i])) ) / s_b

        return hH113, F111


    def calc_induced_source_potential(self, P):
        """Calculates the source-induced potential at the given point P.
        
        Parameters
        ----------
        P : ndarray
            Evaluation point.
        """

        # Calculate geometry
        h,_,_,_,_,a,_,_ = self._calc_geometry(P)

        # Calculate integrals
        hH113, F111 = self._calc_integrals(P)

        phi_s = 0.5*self.sigma*(sum(a*F111) + h*sum(hH113))/np.pi

        return phi_s, hH113, F111, a


    def calc_induced_doublet_potential(self, P):
        """Calculates the doublet-induced potential at the given point P.
        
        Parameters
        ----------
        P : ndarray
            Evaluation point.
        """

        # Calculate integrals
        hH113,_ = self._calc_integrals(P)

        phi_d = 0.5*self.mu*sum(hH113)/np.pi

        return phi_d


class SuperinclinedPanel:
    """A class defining a superinclined constant-doublet-constant-source panel of unit strength in supersonic flow (M_0=sqrt(2)).
    The flow is aligned with the x axis, and the panel lies in the y-z plane.
    
    Parameters
    ----------
    verts : ndarray
        2xN array, the columns of which are the corner points of the panel, expressed in local coords.

    mu : float
        Doublet strength.

    sigma : float
        Source strength.
    """

    def __init__(self, verts, mu, sigma):

        # Store
        self.verts = verts
        self.N = self.verts.shape[1]
        self.mu = mu
        self.sigma = sigma

        # Initialize storage
        self.dx = np.zeros(self.N)
        self.dy = np.zeros(self.N)
        self.t_hat = np.zeros((2,self.N))
        self.n_hat = np.zeros((2,self.N))

        # Calculate edge parameters
        for i in range(self.N):

            # Get index of end vertex
            i_next = (i+1)%self.N

            # Get edge displacements
            self.dx[i] = verts[0,i_next] - verts[0,i]
            self.dy[i] = verts[1,i_next] - verts[1,i]

            # Get edge tangent vector
            self.t_hat[0,i] = self.dx[i]
            self.t_hat[1,i] = self.dy[i]
            self.t_hat[:,i] /= np.linalg.norm(self.t_hat[:,i])

        # Get edge normal vector
        self.n_hat[0,:] = self.t_hat[1,:]
        self.n_hat[1,:] = -self.t_hat[0,:]

        # Calculate area
        if self.N == 3:
            self.A = 0.5*abs(self.dx[0]*self.dy[1] - self.dy[0]*self.dx[1])
        else:
            self.A = abs(self.dx[0]*self.dy[1] - self.dy[0]*self.dx[1])

        # Calculate centroid
        self.c = np.sum(self.verts, axis=1)/self.N


    def _calc_geometry(self, P):
        # Calculates necessary geometry

        # Initialize
        h = P[0]
        a = np.zeros(self.N)
        l1 = np.zeros(self.N)
        l2 = np.zeros(self.N)
        in_dod = np.zeros(self.N, dtype=bool)
        corner_in_dod = np.zeros(self.N, dtype=bool)
        in_dod[:] = False
        corner_in_dod[:] = False

        # Loop through edges
        for i in range(self.N):

            # Get displacement to start vertex
            d = self.verts[:,i] - P[1:]

            # Check if the start vertex is in the DoD
            if d[0]**2 + d[1]**2 < h**2:
                corner_in_dod[i] = True

            # Calculate a
            a[i] = inner2(d, self.n_hat[:,i])
            
            # Calculate l1
            l1[i] = inner2(d, self.t_hat[:,i])

            # Get displacement to end vertex
            i_next = (i+1)%self.N
            d = self.verts[:,i_next] - P[1:]

            # Calculate l2
            l2[i] = inner2(d, self.t_hat[:,i])

        # Loop through edges to check DoD
        for i in range(self.N):

            # Get index of end vertex
            i_next = (i+1)%self.N

            # Check for at least one endpoint in
            if corner_in_dod[i] or corner_in_dod[i_next]:
                in_dod[i] = True
            elif abs(a[i]) < h and h != 0.0 and l1[i]*l2[i] <= 0.0:
                in_dod[i] = True

        return h, a, l1, l2, in_dod


    def _calc_integrals(self, P):
        # Calculates the needed integrals

        # Calculate geometry
        h, a, l1, l2, in_dod = self._calc_geometry(P)


    def calc_induced_source_potential(self, P):
        """Calculates the source-induced potential at the given point P.
        
        Parameters
        ----------
        P : ndarray
            Evaluation point.
        """

        # Get integrals
        self._calc_integrals(P)

        return 0.0, np.zeros(self.N), np.zeros(self.N), np.zeros(self.N)


    def calc_induced_doublet_potential(self, P):
        """Calculates the doublet-induced potential at the given point P.
        
        Parameters
        ----------
        P : ndarray
            Evaluation point.
        """

        # Get integrals
        self._calc_integrals(P)
        
        return 0.0