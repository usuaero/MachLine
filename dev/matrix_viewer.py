import sys
import numpy as np
import matplotlib.pyplot as plt


if __name__=="__main__":

    # Get matrix file
    filename = sys.argv[-1]
    
    # Read into array
    elements = np.genfromtxt(filename)
    element_mags = np.log10(np.abs(elements))
    np.nan_to_num(element_mags, copy=False, neginf=-20.0)
    print(element_mags.shape)

    # Plot
    plt.figure()
    plt.imshow(element_mags, cmap='hot')
    plt.colorbar()
    plt.show()