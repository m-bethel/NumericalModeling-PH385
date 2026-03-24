import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

def plot_comparison():
    # Define the specific files for n=2
    shoot_file = "shoot_n2.csv"
    match_file = "match_n2.csv"
    
    # Check if files exist before proceeding
    if not os.path.exists(shoot_file) or not os.path.exists(match_file):
        print("Missing n=2 data files. Please run your C++ code first.")
        return

    # Create a side-by-side comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    def process_and_plot(filename, ax, title):
        # Load the data: columns are [state, x, psi]
        df = pd.read_csv(filename, names=['state', 'x', 'psi'])
        
        # Get unique states (0 through 10)
        states = sorted(df['state'].unique())
        
        for s in states:
            state_data = df[df['state'] == s]
            x = state_data['x'].values
            psi = state_data['psi'].values

            # Normalize the wavefunction so the total probability is 1
            norm = np.sqrt(np.trapz(psi**2, x)*2)
            if norm > 1E-15:
                psi = psi / norm

            # Parity Mirroring Logic
            # Even states (0, 2, 4...) are symmetric; Odd states (1, 3, 5...) are anti-symmetric
            parity = 1 if int(s) % 2 == 0 else -1
            
            x_full = np.concatenate([-x[::-1], x])
            psi_full = np.concatenate([parity * psi[::-1], psi])
            
            ax.plot(x_full, psi_full, label=f"State {int(s)}")

        ax.set_title(f"{title} (n=2)")
        ax.set_xlabel("Position (x)")
        ax.axhline(0, color='black', lw=0.5)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.85,0.85)
        #ax.set_xlim(-7.5,7.5)

    # Execute plotting for both methods
    process_and_plot(shoot_file, ax1, "Shooting Method")
    process_and_plot(match_file, ax2, "Matching Method")

    ax1.set_ylabel("Wavefunction (psi)")
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="States")
    
    plt.tight_layout()
    plt.savefig("n2_comparison_plot.png")
    print("Saved: n2_comparison_plot.png")
    plt.show()

if __name__ == "__main__":
    plot_comparison()