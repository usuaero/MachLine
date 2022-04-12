import numpy as np

class Source:
    """Defines a supersonic point source.
    
    Parameters
    ----------
    loc : ndarray
        Position.

    sigma : float
        Strength.
    """


    def __init__(self, loc, sigma):

        # Store
        self.loc = loc
        self.sigma = sigma


    def calc_induced_potential(self, P):
        """Calculates the potential induced by this source at point P.
        
        Parameters
        ----------
        P : ndarray
            Influenced point.
        """

        # Check downstream
        phi_s = 0.0
        if P[0] > self.loc[0]:

            # Check Mach cone
            RB2 = (P[0]-self.loc[0])**2 - (P[1]-self.loc[1])**2 - (P[2]-self.loc[2])**2
            if RB2 > 0:

                phi_s = self.sigma/(2.0*np.pi*np.sqrt(RB2))

        return phi_s


class Doublet:
    """Defines a supersonic point doublet aligned with the positive z-axis.
    
    Parameters
    ----------
    loc : ndarray
        Position.
    
    mu : float
        Strength.
    """


    def __init__(self, loc, mu):

        # Store
        self.loc = loc
        self.mu = mu


    def calc_induced_potential(self, P):
        """Calculates the potential induced by this doublet at point P.
        
        Parameters
        ----------
        P : ndarray
            Influenced point.
        """

        # Check downstream
        phi_d = 0.0
        if P[0] > self.loc[0]:

            # Check Mach cone
            RB2 = (P[0]-self.loc[0])**2 - (P[1]-self.loc[1])**2 - (P[2]-self.loc[2])**2
            if RB2 > 0:

                phi_d = -self.mu*(self.loc[2]-P[2])/(2.0*np.pi*RB2**1.5)

        return phi_d