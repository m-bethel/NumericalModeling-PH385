import numpy as np                      # Numerical library for array math
import matplotlib.pyplot as plt         # Plotting library
from matplotlib.animation import FFMpegWriter # The engine that encodes the MP4 video
import glob                             # For finding files using wildcards
import os                               # For file and path operations

# ==========================================
# 1. CONFIGURATION (Must match C++ exactly)
# ==========================================
nx, ny = 1024, 256                      # Grid resolution (X and Y)
Lx, Ly = 4.0, 1.0                       # Physical dimensions of the domain

# --- RENDERING MODE ---
# Set to "speed" for Velocity Magnitude or "pressure" for Pressure
MODE = "pressure" 

# Locate all simulation output files
files = sorted(glob.glob("frame_*.bin"), key=lambda x: int(x.split('_')[1].split('.')[0]))
num_files = len(files)                  # Count how many frames we have to process

if num_files == 0:                      # Safety check if no data is found
    print("Error: No binary files found. Run the C++ simulation first.")
    exit()

# ==========================================
# 2. PLOT & WRITER SETUP
# ==========================================
fig, ax = plt.subplots(figsize=(16, 4), dpi=100) # Create a wide figure at 100 DPI

# Choose colormap: 'magma' is great for speed, 'jet' for pressure
cmap_choice = 'magma' if MODE == "speed" else 'jet'

# Initialize the image object with zeros
im = ax.imshow(np.zeros((ny, nx)), origin='lower', extent=[0, Lx, 0, Ly], 
               cmap=cmap_choice, aspect='equal')

plt.colorbar(im, label=MODE.capitalize()) # Add the color scale bar
ax.set_title(f"CFD Simulation: {MODE.capitalize()}") # Set top title
ax.set_facecolor('black')               # Set background to black for the cylinder hole

# Pre-calculate the cylinder mask (Center at 0.3, 0.49 to match C++ symmetry break)
x_coords = np.linspace(0, Lx, nx)       # X-axis points
y_coords = np.linspace(0, Ly, ny)       # Y-axis points
X, Y = np.meshgrid(x_coords, y_coords)  # Create 2D coordinate grid
mask = (X - 0.3)**2 + (Y - 0.49)**2 <= 0.1**2 # Logical mask for the cylinder

# Configure the MP4 Writer (High bitrate ensures no "blocks" in the gradients)
writer = FFMpegWriter(fps=30, metadata=dict(artist='Matplotlib'), bitrate=8000)

# ==========================================
# 3. THE RENDERING LOOP
# ==========================================
video_filename = f"simulation_{MODE}.mp4" # Define output filename
print(f"Starting render of {num_files} frames to {video_filename}...")

with writer.saving(fig, video_filename, dpi=100): # Open the video file for writing
    for i, filename in enumerate(files):          # Loop through every found binary file
        
        # Read data as 32-bit floats (matches your updated C++ float switch)
        data = np.fromfile(filename, dtype=np.float32)
        
        # Skip files that aren't the right size (prevents crashing if simulation is still running)
        if data.size != 3 * nx * ny:
            continue
            
        # Slicing the 1D data into 2D chunks: [U, V, P]
        u_flat = data[0 : nx*ny]        # First chunk is X-velocity
        v_flat = data[nx*ny : 2*nx*ny]  # Second chunk is Y-velocity
        p_flat = data[2*nx*ny : 3*nx*ny]# Third chunk is Pressure
        
        if MODE == "speed":
            # Physics: Speed is the magnitude of the velocity vector: $\sqrt{u^2 + v^2}$
            grid = np.sqrt(u_flat**2 + v_flat**2).reshape((ny, nx))
            v_min, v_max = 0.0, 2.0     # Typical velocity range for this setup
        else:
            # Just use the pressure grid
            grid = p_flat.reshape((ny, nx))
            v_min, v_max = -0.5, 0.2    # Typical pressure range
            
        grid[mask] = np.nan             # Set the cylinder area to NaN (transparent)
        
        im.set_array(grid)              # Update the image pixels
        im.set_clim(v_min, v_max)       # Apply fixed color limits for consistency
        
        writer.grab_frame()             # Take a "snapshot" of the figure and add to MP4
        
        if i % 20 == 0:                 # Every 20 frames, print progress to console
            print(f"Render Progress: {int(100 * i / num_files)}%")

print(f"Render Complete! Output saved as {video_filename}")