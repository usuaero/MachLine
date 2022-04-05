import numpy as np

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

        # Calculate edge slopes
        self.m = (np.roll(verts[1,:], 1)-verts[1,:])/(np.roll(verts[0,:], 1)-verts[0,:])
        self.l = 1./self.m

        # Calculate area
        d0 = self.verts[:,2] - self.verts[:,1]
        d1 = self.verts[:,1] - self.verts[:,0]
        if self.N == 3:
            self.A = 0.5*abs(d0[0]*d1[1] - d0[1]*d1[0])
        else:
            self.A = abs(d0[0]*d1[1] - d0[1]*d1[0])

        # Calculate centroid
        self.c = np.sum(self.verts, axis=1)/self.N


    def _calc_geometry(self, P):
        # Calculates necessary geometry

        # Set up geometry
        h = P[2]
        xm = np.zeros(self.N)
        ym1 = np.zeros(self.N)
        ym2 = np.zeros(self.N)
        R1 = np.zeros(self.N)
        R2 = np.zeros(self.N)

        # Loop through edges
        for i in range(self.N):
            i_next = (i+1)%self.N

            # Get displacements
            sm1 = P[0] - self.verts[0,i]
            s1 = P[1] - self.verts[1,i]
            sm2 = P[0] - self.verts[0,i_next]
            s2 = P[1] - self.verts[1,i_next]

            # Get hyperbolic radii
            x = sm1**2 - s1**2 - h**2
            if x > 0.0:
                R1[i] = np.sqrt(x)
            else:
                R1[i] = 0.0

            x = sm2**2 - s2**2 - h**2
            if x > 0.0:
                R2[i] = np.sqrt(x)
            else:
                R2[i] = 0.0

            # Calculate edge-based coordinates

            # Subinclined edge
            if abs(self.m[i]) <= 1.0:

                xm[i] = s1 - sm1*self.m[i]
                ym1[i] = sm1 - s1*self.m[i]
                ym2[i] = sm2 - s2*self.m[i]

            # Superinclined edge
            else:

                xm[i] = sm1 - s1*self.l[i]
                ym1[i] = s1 - sm1*self.l[i]
                ym2[i] = s2 - sm2*self.l[i]

        return h, xm, ym1, ym2, R1, R2


    def _calc_integrals(self, P):
        # Calculates the needed integrals

        # Calculate geometry
        h, xm, ym1, ym2, R1, R2 = self._calc_geometry(P)

        # Initialize
        Q1 = 0.0
        w0 = np.zeros(self.N)

        # Loop through edges
        for i in range(self.N):

            # Calculate Q1

            # Subsonic or sonic edge
            if abs(self.m[i]) <= 1.0:

                # Both endpoints in
                if R1[i] != 0.0 and R2[i] != 0.0:

                    A = h*xm[i]*(ym1[i]*R2[i] - ym2[i]*R1[i])
                    B = h**2*ym1[i]*ym2[i] + xm[i]**2*R1[i]*R2[i]

                    Q1 += np.arctan2(A, B)

                # First endpoint in
                elif R1[i] != 0.0:

                    Q1 += -np.sign(h)*np.arctan2(xm[i]*R1[i], abs(h)*ym1[i])

                # Second endpoint in
                elif R1[i] != 0.0:

                    Q1 += np.sign(h)*np.arctan2(xm[i]*R2[i], abs(h)*ym2[i])

            # Calculate w0
            
            # Subsonic edge
            if abs(self.m[i]) < 1.0:

                x = np.sqrt(1.0-self.m[i]**2)

                # Both endpoints in
                if R1[i] != 0.0 and R2[i] != 0.0:

                    A = ym2[i] + x*R2[i]
                    B = ym1[i] + x*R1[i]

                    w0[i] = 1.0/x*np.log(A/B)

        return Q1, w0


    def calc_induced_source_potential(self, P):
        """Calculates the source-induced potential at the given point P.
        
        Parameters
        ----------
        P : ndarray
            Evaluation point.
        """

        # Calculate geometry
        h, xm, ym1, ym2, R1, R2 = self._calc_geometry(P)

        # Calculate integrals
        Q1, w0 = self._calc_integrals(P)

        # Loop through edges to add up w0
        phi_s = 0.0
        for i in range(self.N):

            if abs(self.m[i]) < 1.0:
                phi_s += xm[i]*w0[i]*self.m[i]
            else:
                phi_s += xm[i]*w0[i]

        # Add in Q1
        phi_s = 0.5*self.sigma*(phi_s - h*Q1)/np.pi

        return phi_s


    def calc_induced_doublet_potential(self, P):
        """Calculates the doublet-induced potential at the given point P.
        
        Parameters
        ----------
        P : ndarray
            Evaluation point.
        """

        # Calculate geometry
        h, xm, ym1, ym2, R1, R2 = self._calc_geometry(P)

        # Calculate integrals
        Q1, w0 = self._calc_integrals(P)

        phi_d = -0.5*self.mu*Q1/np.pi

        return phi_d