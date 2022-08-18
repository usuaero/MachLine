import sys
import numpy as np
import matplotlib.pyplot as plt


if __name__=="__main__":

    # Get matrix file
    filename = sys.argv[-1]
    
    # Read into array
    elements = np.genfromtxt(filename, dtype=np.single)
    elements = np.abs(elements, dtype=np.single)
    elements = np.log10(elements, dtype=np.single)
    np.nan_to_num(elements, copy=False, neginf=-20.0)
    print(elements.shape)

    # Plot
    plt.figure()
    plt.imshow(elements, cmap='hot')
    plt.colorbar()
    plt.show()