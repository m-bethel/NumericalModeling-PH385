import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm # Colormap library for the gradient styling

# 1. Determine which specific data file to load based on the omega you want to visualize
target_omega = 1.9
# Formats the string to match the C++ std::to_string output (which usually adds 6 decimal places)
filename = f"potential_w_{target_omega:.6f}.dat" 

# Load the raw text file into a flat numpy array structure
data = np.loadtxt(filename)

# Slice the data into three separate 1D arrays for X, Y, and Voltage
x = data[:, 0]
y = data[:, 1]
v = data[:, 2]

# 2. Mathematically determine the grid resolution (e.g., sqrt of 1681 points = 41)
grid_size = int(np.sqrt(len(x)))
print(grid_size)

# Reshape the flat 1D arrays into 2D matrices required by Matplotlib's 3D surface plotter
X = x.reshape((grid_size, grid_size))
Y = y.reshape((grid_size, grid_size))
V = v.reshape((grid_size, grid_size))

# 3. Initialize the plotting window and define the 3D axes
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d') # '111' means 1x1 grid, 1st subplot

# Plot the surface. 
# cmap=cm.viridis applies the blue/yellow color scale.
# antialiased=True smooths the visual artifacts on the mesh lines.
surf = ax.plot_surface(X, Y, V, cmap=cm.viridis,
                       linewidth=0.1, antialiased=True, edgecolors='black', alpha=0.8)

# 4. Lock the vertical (Z) axis to exactly -0.5 to 0.5. 
# This prevents the camera from auto-scaling and makes it easier to compare different graphs.
ax.set_zlim(-0.5, 0.5) 

# Add aesthetic labels
ax.set_title('Electric Potential at z=0')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('V (Volts)')

# Draw the legend bar on the right side showing which colors map to which voltage values
fig.colorbar(surf, shrink=0.5, aspect=5)

# Render the window to the screen
plt.show()