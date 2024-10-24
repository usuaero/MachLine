import pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

# # Set global font properties for the entire figure
# plt.rcParams['font.family'] = 'serif'
# plt.rcParams['font.serif'] = ['Times New Roman']
# plt.rcParams['font.size'] = 14

# Load figure from a pickled file (.pkl)
with open("C:\\Users\\nathan\\git-repos\\MachLine\\studies\\adjoint_studies\\super_11\pickle_jar\\super_octa_5_cp_offsets_dCFx.pkl", 'rb') as f:
    fig = pickle.load(f)  # fig is now a Matplotlib Figure object

num_steps = 5

# Access the axes
ax = fig.gca()  # Get current axes

# Iterate over each line in the plot
for line in ax.get_lines():
    # Get the data (x and y values)
    x_data = line.get_xdata()
    y_data = line.get_ydata()

    # Remove the first data point
    line.set_xdata(x_data[1:])  # Remove the first element from x data
    line.set_ydata(y_data[1:])  # Remove the first element from y data



max_dash = 10  # Max dash length
min_dash = 1   # Min dash length
dash_step = (max_dash - min_dash) / (num_steps - 1)
# Create dashes starting from max_dash to min_dash
dash_styles = [(max_dash - i * dash_step, 2) for i in range(num_steps)]
return dash_styles

    # Add markers (shapes) at each data point
    # ax.scatter(x_data[1:], y_data[1:], marker='o', color='k', s=10)  # Adjust marker type, color, and size as needed

# Set the y-axis label with F and subscript x
ax.set_ylabel("Norm of CFx Vector", fontsize=14, labelpad=10)  # F followed by Unicode subscript x

ax.set_xlabel("Control Point Offset", fontsize=14, labelpad=7)  # F followed by Unicode subscript x

# Set the title with F and subscript x
ax.set_title("Norm of CFx vs Control Point Offset", fontsize=14)  # F followed by Unicode subscript x

#
# new_x_data9 = np.array([1e-10, 1e-9, 1e-8, 1e-7])
# new_y_data9 = np.array([0.1, 0.05, 0.03, 0.01])  # Replace with your actual data

# # Add the new data series to the plot with superscript in the label
# ax.plot(new_x_data9, new_y_data9, label="Step size = $10^{-5}$", color='blue', marker='o')


# Set x-axis to logarithmic scale
ax.set_xscale('log')
ax.set_xlim(left=5.0e-11)  # Set the left limit to 1.5 * 10^-10

# Set the number of ticks on the logarithmic scale
ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=10))

# Display the modified figure
plt.show()

# Optional: Save the modified figure back to a pickle file
with open("C:\\Users\\nathan\\git-repos\\MachLine\\studies\\adjoint_studies\\central_diff_norm_cp_offset\\pickle_jar\\modified_figure_with_shapes.pkl", 'wb') as f:
    pickle.dump(fig, f)