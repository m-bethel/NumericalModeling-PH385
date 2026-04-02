import matplotlib.pyplot as plt
import numpy as np

# Load the data from your older C++ output method
# Assuming the file has x, y, V columns separated by whitespace/tabs

# 1. Determine which specific data file to load based on the omega you want to visualize
target_omega = 1.9
# Formats the string to match the C++ std::to_string output (which usually adds 6 decimal places)
filename = f"potential_w_{target_omega:.6f}.dat" 

# Load the raw text file into a flat numpy array structure
data = np.loadtxt(filename)

# Extract columns into 1D arrays
x = data[:, 0]
y = data[:, 1]
v = data[:, 2]

# Calculate grid dimension (41 for your current 0.025 spacing setup)
grid_size = int(np.sqrt(len(x)))

# Reshape the 1D arrays into the 2D matrices necessary for contour mapping
X = x.reshape((grid_size, grid_size))
Y = y.reshape((grid_size, grid_size))
V = v.reshape((grid_size, grid_size))

# Create the figure window
plt.figure(figsize=(8, 6))

# Generate filled contour plots. 
# 'levels=50' tells it to draw 50 distinct color bands to make the gradient look smooth.
cp = plt.contourf(X, Y, V, levels=50, cmap='viridis')

# Attach the color legend to the side of the plot
plt.colorbar(cp, label='Potential (V)')

# Add standard labels to the axes
plt.title('Electric Potential at z=0')
plt.xlabel('x (m)')
plt.ylabel('y (m)')

# Render the window to the screen
plt.show()