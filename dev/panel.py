import numpy as np

class Panel:
    """A class defining a constant-doublet-constant-source panel of unit strength in supersonic flow (M_0=sqrt(2)).
    
    Parameters
    ----------
    verts : ndarray
        2x3 array, the columns of which are the corner points of the panel, expressed in local coords.
    """

    def __init__(self, verts):

        # Store
        self.verts = verts

        # Calculate edge slopes
        self.m = (np.roll(verts[1,:], 1)-verts[1,:])/(np.roll(verts[0,:], 1)-verts[0,:])
        self.l = 1./self.m

        # Calculate area
        d0 = self.verts[:,2] - self.verts[:,1]
        d1 = self.verts[:,1] - self.verts[:,0]
        self.A = 0.5*np.abs(d0[0]*d1[1] - d0[1]*d1[0])