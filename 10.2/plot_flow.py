import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

fig, ax = plt.subplots()

def update(frame):
    ax.clear()
    data = pd.read_csv(f'frame_{frame}.csv')

    num_points = len(data)
    #print(f"Points in CSV: {num_points}")

    nx, ny = 32,16

    if num_points != nx * ny:
        print(f"ERROR: CSV has {num_points} points, but we expected {nx*ny}!")
        return

    X = data['x'].values.reshape(ny, nx)
    Y = data['y'].values.reshape(ny, nx)
    P = data['p'].values.reshape(ny, nx)
    U = data['u'].values.reshape(ny, nx)
    V = data['v'].values.reshape(ny, nx)

    cont = ax.contourf(X, Y, P, cmap='RdBu_r', levels=np.linspace(np.min(P), np.max(P), 20))
    #cont = ax.contourf(X, Y, speed, cmap='jet', levels=50, vmax=0.3)
    #ax.quiver(X[::2, ::2], Y[::2, ::2], U[::2, ::2], V[::2, ::2], color='black', alpha=0.5)
    
    ax.set_title(f"Lid-Driven Cavity: Frame {frame}")
    return cont

#animation = ani.FuncAnimation(fig, update, frames=range(0, 150), interval=50)
animation1 = ani.FuncAnimation(fig, update, frames=range(0, 3), interval=500)
plt.show()