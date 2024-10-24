import pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker


# Load figure from a pickled file (.pkl)
with open("C:\\Users\\nathan\\git-repos\\MachLine\\studies\\adjoint_studies\\super_11\pickle_jar\\super_octa_5_cp_offsets_dCFx.pkl", 'rb') as f:
    fig = pickle.load(f)  # fig is now a Matplotlib Figure object

num_steps = 5

# Access the axes
ax = fig.gca()  # Get current axes

plt.show()