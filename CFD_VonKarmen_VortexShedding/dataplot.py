import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import os

# ==========================================
# 1. CONFIGURATION
# ==========================================
nx, ny = 1024, 256 
Lx, Ly = 4.0, 1.0

# --- VISUALIZATION MODE ---
# Set this to "speed" to see Velocity Magnitude, or "pressure" for Pressure
MODE = "speed" 

num_files = len(glob.glob("frame_*.bin"))
if num_files == 0:
    print("Error: No binary files found.")
    exit()

# ==========================================
# 2. PLOT SETUP
# ==========================================
fig, ax = plt.subplots(figsize=(16, 4))

# Use 'magma' or 'inferno' for Speed (dark to bright)
# Use 'jet' or 'RdBu_r' for Pressure
cmap_choice = 'magma' if MODE == "speed" else 'terrain'

im = ax.imshow(np.zeros((ny, nx)), origin='lower', extent=[0, Lx, 0, Ly], 
               cmap=cmap_choice, animated=True, aspect='equal')

plt.colorbar(im, label=MODE.capitalize())
ax.set_title(f"Fluid {MODE.capitalize()}")
ax.set_facecolor('black')

# Pre-calculate coordinate grid for the cylinder mask
x_coords = np.linspace(0, Lx, nx)
y_coords = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x_coords, y_coords)
# Cylinder offset at 0.49 to match your C++ symmetry break
mask = (X - 0.3)**2 + (Y - 0.49)**2 <= 0.1**2

# ==========================================
# 3. UPDATE FUNCTION (WITH SPEED LOGIC)
# ==========================================
def update(frame):
    filename = f"frame_{frame}.bin"
    if not os.path.exists(filename): return [im]
        
    # Read binary as float32 (matches your C++ 'float' change)
    data = np.fromfile(filename, dtype=np.float32)
    
    if data.size != 3 * nx * ny: return [im]
    
    # --- DATA SLICING ---
    # The file contains U, then V, then P. Each is (nx * ny) long.
    u_flat = data[0 : nx*ny]
    v_flat = data[nx*ny : 2*nx*ny]
    p_flat = data[2*nx*ny : 3*nx*ny]
    
    if MODE == "speed":
        # Calculate Magnitude: sqrt(u^2 + v^2)
        # This shows where the fluid is moving fastest vs slowest
        grid = np.sqrt(u_flat**2 + v_flat**2).reshape((ny, nx))
        v_min, v_max = 0.0, 2.0  # Speed range (Inlet is 1.45)
    else:
        # Just show Pressure
        grid = p_flat.reshape((ny, nx))
        v_min, v_max = -0.5, 0.5 # Typical pressure range
    
    # Apply the black cylinder mask
    grid[mask] = np.nan
    
    # Update image data
    im.set_array(grid)
    im.set_clim(v_min, v_max)
    
    return [im]

# ==========================================
# 4. RUN
# ==========================================
ani = FuncAnimation(fig, update, frames=range(0, num_files), 
                    blit=True, interval=1, repeat=True, cache_frame_data=False)

plt.tight_layout()
plt.show()