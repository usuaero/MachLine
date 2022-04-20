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
                l1[i] = -np.sqrt(abs(g2[i])) # From PAN AIR
                R1[i] = 0.0

            x = (g2[i] - l2[i]**2)/self.b[i]
            if x > 0.0 and d2[0] < 0.0:
                R2[i] = np.sqrt(x)
            else:
                l2[i] = np.sqrt(abs(g2[i])) # From PAN AIR
                R2[i] = 0.0

            # Check DoD

            # If the point is upstream of the edge, it is not in
            if d1[0] > 0.0 and d2[0] > 0.0:
                in_dod[i] = False

            # Check if both endpoints are out
            elif R1[i] == 0.0 and R2[i] == 0.0:

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


        return h, R1, R2, l1, l2, a, g2, in_dod


    def _calc_integrals(self, P):
        # Calculates the needed integrals

        # Calculate geometry
        h, R1, R2, l1, l2, a, g2, in_dod = self._calc_geometry(P)

        # Initialize
        hH113 = 0.0
        F111 = np.zeros(self.N)
        print()
        print(P)
        print(in_dod)

        # Loop through edges
        for i in range(self.N):

            # Check DoD
            if in_dod[i]:

                # Calculate hH(1,1,3)
                s_b = np.sqrt(abs(self.b[i]))

                # Check for point on panel plane
                if h != 0.0:

                    # Supersonic edge
                    if self.b[i] > 0.0:

                        # Neither endpoint in
                        if R1[i] == 0.0 and R2[i] == 0.0:

                            # Calculate hH113
                            hH113 += np.pi*np.sign(h*self.n_hat[0,i])

                        # At least one endpoint in
                        else:
                            F1 = (l1[i]*R2[i] - l2[i]*R1[i]) / g2[i]
                            F2 = (self.b[i]*R1[i]*R2[i] + l1[i]*l2[i]) / g2[i]
                            hH113 += np.arctan2(h*a[i]*F1, R1[i]*R2[i] + h**2*F2)
                
                    # Subsonic edge
                    else:

                        # Caculate preliminaries
                        F1 = (R2[i]**2 - R1[i]**2) / (l1[i]*R2[i] + l2[i]*R1[i])
                        F2 = (g2[i] - l1[i]**2 - l2[i]**2) / (self.b[i]*R1[i]*R2[i] - l1[i]*l2[i])

                        # Add for edge
                        hH113 += np.arctan2(h*a[i]*F1, R1[i]*R2[i] + h**2*F2)


                # Calculate F(1,1,1)

                # Supersonic edge
                if self.b[i] > 0.0:

                    # Check for Mach wedge condition
                    if R1[i] == 0.0 and R2[i] == 0.0:

                        F111[i] = np.pi/s_b

                    # At least one in
                    else:

                        F1 = (l1[i]*R2[i] - l2[i]*R1[i]) / g2[i]
                        F2 = (self.b[i]*R1[i]*R2[i] + l1[i]*l2[i]) / g2[i]
                        F111[i] = -np.arctan2(s_b*F1, F2) / s_b

                # Subsonic edge
                else:

                    F1 = s_b*R1[i] + abs(l1[i])
                    F2 = s_b*R2[i] + abs(l2[i])
                    F111[i] = -np.sign(self.n_hat[1,i]) * np.log(F1/F2) / s_b

        # Check for on panel
        if h == 0.0 and all(a > 0.0):
            hH113 = 2.0*np.pi

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

        phi_s = 0.5*self.sigma*(sum(a*F111) + h*hH113)/np.pi

        return phi_s, hH113, F111


    def calc_induced_doublet_potential(self, P):
        """Calculates the doublet-induced potential at the given point P.
        
        Parameters
        ----------
        P : ndarray
            Evaluation point.
        """

        # Calculate integrals
        hH113,_ = self._calc_integrals(P)

        phi_d = 0.5*self.mu*hH113/np.pi

        return phi_d


#class Panel:
#    """A class defining a constant-doublet-constant-source panel of unit strength in supersonic flow (M_0=sqrt(2)).
#    
#    Parameters
#    ----------
#    verts : ndarray
#        2xN array, the columns of which are the corner points of the panel, expressed in local coords.
#
#    mu : float
#        Doublet strength.
#
#    sigma : float
#        Source strength.
#    """
#
#    def __init__(self, verts, mu, sigma):
#
#        # Store
#        self.verts = verts
#        self.N = self.verts.shape[1]
#        self.mu = mu
#        self.sigma = sigma
#
#        # Initialize storage
#        self.m = np.zeros(self.N)
#        self.l = np.zeros(self.N)
#        self.dx = np.zeros(self.N)
#        self.dy = np.zeros(self.N)
#        self.t_hat = np.zeros((2,self.N))
#        self.n_hat = np.zeros((2,self.N))
#        self.b = np.zeros(self.N)
#
#        # Calculate edge parameters
#        for i in range(self.N):
#
#            # Get index of end vertex
#            i_next = (i+1)%self.N
#
#            # Get edge displacements
#            self.dx[i] = verts[0,i_next] - verts[0,i]
#            self.dy[i] = verts[1,i_next] - verts[1,i]
#
#            # Get edge slopes
#            self.m[i] = self.dy[i]/self.dx[i]
#
#            # Get edge tangent vector
#            self.t_hat[0,i] = self.dx[i]
#            self.t_hat[1,i] = self.dy[i]
#            self.t_hat[:,i] /= np.linalg.norm(self.t_hat[:,i])
#
#        # Get edge normal vector
#        self.n_hat[0,:] = self.t_hat[1,:]
#        self.n_hat[1,:] = -self.t_hat[0,:]
#
#        # Get b
#        self.b = self.n_hat[0,:]**2 - self.n_hat[1,:]**2
#
#        # Calculate edge inverse slope
#        self.l = 1./self.m
#
#        # Calculate area
#        if self.N == 3:
#            self.A = 0.5*abs(self.dx[0]*self.dy[1] - self.dy[0]*self.dx[1])
#        else:
#            self.A = abs(self.dx[0]*self.dy[1] - self.dy[0]*self.dx[1])
#
#        # Calculate centroid
#        self.c = np.sum(self.verts, axis=1)/self.N
#
#
#    def _calc_geometry(self, P):
#        # Calculates necessary geometry
#
#        # Set up geometry
#        h = P[2]
#        xm = np.zeros(self.N)
#        ym1 = np.zeros(self.N)
#        ym2 = np.zeros(self.N)
#        R1 = np.zeros(self.N)
#        R2 = np.zeros(self.N)
#        l1 = np.zeros(self.N)
#        l2 = np.zeros(self.N)
#        a = np.zeros(self.N)
#        g2 = np.zeros(self.N)
#
#        # Loop through edges
#        for i in range(self.N):
#            i_next = (i+1)%self.N
#
#            # Get displacements
#            # First vertex
#            d1 = P[0:2] - self.verts[:,i]
#            sm1 = d1[0]
#            s1 = d1[1]
#            l1[i] = -self.n_hat[1,i]*d1[0] - self.n_hat[0,i]*d1[1]
#
#            # Second vertex
#            d2 = P[0:2] - self.verts[:,i_next]
#            sm2 = d2[0]
#            s2 = d2[1]
#            l2[i] = -self.n_hat[1,i]*d2[0] - self.n_hat[0,i]*d2[1]
#
#            # Perpendicular distances
#            a[i] = inner2(self.n_hat[:,i], -d1)
#            g2[i] = a[i]**2 - self.b[i]*h**2
#
#            # Get hyperbolic radii
#            x = sm1**2 - s1**2 - h**2
#            x = (g2[i] - l1[i]**2)/self.b[i]
#            if x > 0.0:
#                R1[i] = np.sqrt(x)
#            else:
#                R1[i] = 0.0
#
#            x = sm2**2 - s2**2 - h**2
#            x = (g2[i] - l2[i]**2)/self.b[i]
#            if x > 0.0:
#                R2[i] = np.sqrt(x)
#            else:
#                R2[i] = 0.0
#
#            # Calculate edge-based coordinates
#
#            # Subsonic or sonic edge (xhatm and yhatm)
#            if abs(self.m[i]) <= 1.0:
#                
#                # Get coordinate directions
#
#                # e_x^m
#                # According to Davis, this axis points towards the direction of integration.
#                # "The axis itself is defined by a rotation of theta from the panel's eta-direction about
#                # the panel's zeta axis which points out of the page; theta here is the angle of the edge
#                # relative to the freestream."
#                e_xm = np.array([self.m[i], -1.0])
#               
#                # Make sure e_x^m points in the direction of integration
#                # Except, under the above definition, xm is perpendicular to the edge, so this is zero.
#                if inner2(e_xm, [self.dx[i], self.dy[i]]) < 0.0:
#                    e_xm = -e_xm
#
#                # e_y^m
#                # According to Davis, this axis is the axis of the edge itself.
#                # This points toward the negative xi-direction of the panel coordinate system.
#                e_ym = np.array([-1.0, -self.m[i]])
#
#                # Get coordinates
#                xm[i] = inner2(e_xm, d1)
#                ym1[i] = inner2(e_ym, d1)
#                ym2[i] = inner2(e_ym, d2)
#
#                xm_alt = inner2(e_xm, d2)
#                if abs(xm[i] -xm_alt ) > 1e-12:
#                    print("xm not constant for edge!!!!!!!!")
#                    print(xm[i])
#                    print(xm_alt)
#
#            # Supersonic edge
#            else:
#
#                # Get coordinate directions
#                e_xm = np.array([1.0, -self.l[i]])
#                e_ym = np.array([-self.l[i], 1.0])
#
#                # Get coordinates
#                xm[i] = inner2(e_xm, d1)
#                ym1[i] = inner2(e_ym, d1)
#                ym2[i] = inner2(e_ym, d2)
#
#        return h, xm, ym1, ym2, R1, R2, l1, l2, a, g2
#
#
#    def _calc_integrals(self, P):
#        # Calculates the needed integrals
#
#        # Calculate geometry
#        h, xm, ym1, ym2, R1, R2, l1, l2, a, g2 = self._calc_geometry(P)
#
#        # Initialize
#        Q1 = 0.0
#        w0 = np.zeros(self.N)
#        hH113 = 0.0
#        F111 = np.zeros(self.N)
#
#        # Loop through edges
#        for i in range(self.N):
#
#            # Calculate Q1
#
#            # Subsonic or sonic edge
#            if abs(self.m[i]) <= 1.0:
#
#                # Both endpoints in
#                if R1[i] != 0.0 and R2[i] != 0.0:
#
#                    A = h*xm[i]*(ym1[i]*R2[i] - ym2[i]*R1[i])
#                    B = h**2*ym1[i]*ym2[i] + xm[i]**2*R1[i]*R2[i]
#
#                    Q1 += np.arctan2(A, B)
#
#                # First endpoint in
#                elif R1[i] != 0.0:
#
#                    Q1 += -np.sign(h)*np.arctan2(xm[i]*R1[i], abs(h)*ym1[i])
#
#                # Second endpoint in
#                elif R2[i] != 0.0:
#
#                    Q1 += np.sign(h)*np.arctan2(xm[i]*R2[i], abs(h)*ym2[i])
#
#            # Calculate w0
#            
#            # Subsonic edge
#            if abs(self.m[i]) < 1.0:
#
#                x = np.sqrt(1.0-self.m[i]**2)
#
#                # Both endpoints in
#                if R1[i] != 0.0 and R2[i] != 0.0:
#
#                    A = ym2[i] + x*R2[i]
#                    B = ym1[i] + x*R1[i]
#
#                    w0[i] = 1.0/x*np.log(A/B)
#
#            # Calculate hH(1,1,3)
#
#            # Caculate preliminaries
#            if self.b[i] >= 0.0:
#                F1 = (l1[i]*R2[i] - l2[i]*R1[i]) / g2[i]
#                F2 = (self.b[i]*R1[i]*R2[i] + l1[i]*l2[i]) / g2[i]
#            else:
#                F1 = (R2[i]**2 - R1[i]**2) / (l1[i]*R2[i] + l2[i]*R1[i])
#                F2 = (g2[i] - l1[i]**2 - l2[i]**2) / (self.b[i]*R1[i]*R2[i] - l1[i]*l2[i])
#
#            # Add for edge
#            hH113 -= np.arctan2(h*a[i]*F1, R1[i]*R2[i] + h**2*F2)
#
#            # Calculate F(1,1,1)
#            if self.b[i] > 0.0:
#                F111[i] = -1.0/np.sqrt(self.b[i]) * np.arctan2(np.sqrt(self.b[i])*F1, F2)
#            else:
#                F1 = np.sqrt(-self.b[i])*R1[i] + abs(l1[i])
#                F2 = np.sqrt(-self.b[i])*R2[i] + abs(l2[i])
#                F111[i] = -np.sign(self.n_hat[1,i])/np.sqrt(-self.b[i]) * np.log(F1/F2)
#
#        return Q1, w0, hH113, F111
#
#
#    def calc_induced_source_potential(self, P):
#        """Calculates the source-induced potential at the given point P.
#        
#        Parameters
#        ----------
#        P : ndarray
#            Evaluation point.
#        """
#
#        # Calculate geometry
#        h, xm, ym1, ym2, R1, R2, l1, l2, a, g2 = self._calc_geometry(P)
#
#        # Calculate integrals
#        Q1, w0, hH113, F111 = self._calc_integrals(P)
#
#        # Loop through edges to add up w0
#        phi_s = 0.0
#        for i in range(self.N):
#
#            if abs(self.m[i]) < 1.0:
#                phi_s += xm[i]*w0[i]*self.m[i]
#            else:
#                phi_s += xm[i]*w0[i]
#
#        # Add in Q1
#        phi_s = 0.5*self.sigma*(phi_s - h*Q1)/np.pi
#
#        phi_s = 0.5*self.sigma*(sum(a*F111) - h*hH113)/np.pi
#
#        return phi_s
#
#
#    def calc_induced_doublet_potential(self, P):
#        """Calculates the doublet-induced potential at the given point P.
#        
#        Parameters
#        ----------
#        P : ndarray
#            Evaluation point.
#        """
#
#        # Calculate geometry
#        h, xm, ym1, ym2, R1, R2, l1, l2, a, g2 = self._calc_geometry(P)
#
#        # Calculate integrals
#        Q1, w0, hH113, F111 = self._calc_integrals(P)
#
#        phi_d = -0.5*self.mu*Q1/np.pi
#
#        phi_d = -0.5*self.mu*hH113/np.pi
#
#        return phi_d