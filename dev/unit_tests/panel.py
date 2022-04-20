import numpy as np


def inner2(x,y):
    return x[0]*y[0] + x[1]*y[1]


class Panel:
    """A class defining a constant-doublet-constant-source panel of unit strength in supersonic flow (M_0=sqrt(2)).
    
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
            x = (g2[i] - l1[i]**2)/self.b[i]
            if x > 0.0 and d1[0] < 0.0:
                R1[i] = np.sqrt(x)
            else:
                R1[i] = 0.0

            x = (g2[i] - l2[i]**2)/self.b[i]
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

                # Supersonic edge
                if self.b[i] > 0.0:

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