
import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV you created in C++
data = pd.read_csv('convergence_data.csv')

# Filter out the 2.0 value if it makes the scale too hard to read
# (Since it didn't actually converge)
converged_data = data[data['omega'] < 2.0]

plt.figure(figsize=(10, 6))
plt.plot(converged_data['omega'], converged_data['iterations'], 'ro-', markersize=8, label='SOR Convergence')

# Find the minimum
min_iters = converged_data['iterations'].min()
opt_omega = converged_data.loc[converged_data['iterations'] == min_iters, 'omega'].values[0]

plt.title('Convergence Study: Finding Optimal $\omega$')
plt.xlabel('$\omega$ (Over-relaxation Parameter)')
plt.ylabel('Number of Iterations')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()