#include "program.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

// Constructor: Initialize the potential and fixed arrays, and define point charges
PotentialSolver::PotentialSolver() : omega(1.5){
    // Resize the 3D potential array to GRID_SIZE in all three dimensions and fill with 0.0 Volts
    V.resize(GRID_SIZE, vector<vector<double>>(GRID_SIZE, vector<double>(GRID_SIZE, 0.0)));
    
    // Resize the 3D fixed boolean array to match, defaulting to false (no nodes locked yet)
    fixed.resize(GRID_SIZE, vector<vector<bool>>(GRID_SIZE, vector<bool>(GRID_SIZE, false)));
    
    // Instantiate the four point charges and add them to the charges vector
    charges.push_back({0.25, 0.0, 0.0, CHARGE});    // +1 uC on the positive x-axis
    charges.push_back({0.0, 0.25, 0.0, CHARGE});    // +1 uC on the positive y-axis
    charges.push_back({-0.25, 0.0, 0.0, -CHARGE});  // -1 uC on the negative x-axis
    charges.push_back({0.0, -0.25, 0.0, -CHARGE});  // -1 uC on the negative y-axis
}

// Set up Dirichlet boundary conditions and the initial potential distribution
void PotentialSolver::initializeConditions() {
    // 1. Set outer boundary conditions: all six outer faces of the simulation box are grounded (0V)
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            fixed[0][i][j] = true;           // Lock -x face
            fixed[GRID_SIZE-1][i][j] = true; // Lock +x face
            fixed[i][0][j] = true;           // Lock -y face
            fixed[i][GRID_SIZE-1][j] = true; // Lock +y face
            fixed[i][j][0] = true;           // Lock -z face
            fixed[i][j][GRID_SIZE-1] = true; // Lock +z face
        }
    }
    
    // 2. Lock the grid points closest to the exact coordinates of the point charges
    for (const auto& charge : charges) {
        // Map physical (x,y,z) coordinates to grid (i,j,k) integer indices
        int i = (charge.x - GRID_MIN) / GRID_SPACING;
        int j = (charge.y - GRID_MIN) / GRID_SPACING;
        int k = (charge.z - GRID_MIN) / GRID_SPACING;

        // Ensure the calculated index is safely within the bounds of our 3D array
        if (i >= 0 && i < GRID_SIZE && j >= 0 && j < GRID_SIZE && k >= 0 && k < GRID_SIZE) {
            // Because true point charges approach infinity at r=0, we cap the potential to +/- 0.5V
            V[i][j][k] = (charge.q > 0) ? 0.5 : -0.5; 
            fixed[i][j][k] = true; // Lock these cells so the solver doesn't overwrite them
        }
    }

    // 3. Construct the flat, conducting disk in the center of the box
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            // Convert current i, j grid indices back into physical x, y distances
            double x = GRID_MIN + i * GRID_SPACING;
            double y = GRID_MIN + j * GRID_SPACING;
            
            // Calculate radial distance from the center z-axis using Pythagorean theorem
            double r = sqrt(x*x + y*y);
            
            // Calculate the exact middle index to place the disk on the z=0 plane
            int k = GRID_SIZE / 2;
            
            // If the grid point falls inside the geometric radius of the disk...
            if (r <= DISK_RADIUS) {
                V[i][j][k] = DISK_VOLTAGE;  // Set the potential to the required 0.25V
                fixed[i][j][k] = true;      // Lock it so it remains a perfect conductor
            }
        }
    }
    
    // 4. Create an educated initial guess for the free space to speed up convergence
    for (int i = 1; i < GRID_SIZE-1; i++) {           // Loop x (skipping walls)
        for (int j = 1; j < GRID_SIZE-1; j++) {       // Loop y (skipping walls)
            for (int k = 1; k < GRID_SIZE-1; k++) {   // Loop z (skipping walls)
                // Only guess values for points that aren't locked boundaries
                if (!fixed[i][j][k]) {
                    // Find physical coordinates of the current free node
                    double x = GRID_MIN + i * GRID_SPACING;
                    double y = GRID_MIN + j * GRID_SPACING;
                    double z = GRID_MIN + k * GRID_SPACING;
                    
                    // Sum the analytical Coulomb potential from all 4 charges (Superposition principle)
                    for (const auto& charge : charges) {
                        double dx = x - charge.x;
                        double dy = y - charge.y;
                        double dz = z - charge.z;
                        
                        double r = sqrt(dx*dx + dy*dy + dz*dz); // Distance to the charge
                        
                        // Prevent dividing by numbers too close to zero to avoid numerical explosions
                        if (r > GRID_SPACING/2) {
                            V[i][j][k] += K * charge.q / r; // Add V = kq/r to the node's total
                        }
                    }
                }
            }
        }
    }
}

// Perform one full sweep of the Successive Over-Relaxation (SOR) algorithm
double PotentialSolver::iterate() {
    double max_diff = 0.0; // Variable to track the largest change made to any single node
    
    // Loop through the entire 3D inner volume (ignoring the fixed boundary walls at index 0 and GRID_SIZE-1)
    for (int i = 1; i < GRID_SIZE - 1; i++) {
        for (int j = 1; j < GRID_SIZE - 1; j++) {
            for (int k = 1; k < GRID_SIZE - 1; k++) {
                
                // If this node is part of the disk or a charge, skip it—do not update
                if (fixed[i][j][k]) continue;

                double old_v = V[i][j][k]; // Store the current potential before modifying it
                
                // Finite difference approximation of Laplace's equation in 3D:
                // The new potential should theoretically be the average of its 6 nearest neighbors
                double target_v = (V[i+1][j][k] + V[i-1][j][k] + 
                                   V[i][j+1][k] + V[i][j-1][k] + 
                                   V[i][j][k+1] + V[i][j][k-1]) / 6.0;

                // Apply the SOR formula: instead of just stepping to the target_v (which is Jacobi/Gauss-Seidel),
                // we 'overshoot' the target by a factor of omega to accelerate convergence.
                V[i][j][k] = old_v + omega * (target_v - old_v);

                // Calculate how much this specific node changed
                double diff = std::abs(V[i][j][k] - old_v);
                // Keep a record of the absolute largest change in the whole grid this sweep
                if (diff > max_diff) max_diff = diff;
            }
        }
    }
    // Return the largest change so the solve() function knows if we hit tolerance
    return max_diff;
}

// Controls the main loop, running iterate() until the system settles into equilibrium
int PotentialSolver::solve() {
    // Output configuration to terminal
    cout << "Starting over-relaxation solver...\n";
    cout << "Grid size: " << GRID_SIZE << "x" << GRID_SIZE << "x" << GRID_SIZE << endl;
    cout << "Omega: " << omega << endl;
    
    int iteration = 0;      
    double max_change;      
    
    // Main execution loop
    do {
        // Run one sweep of the grid and capture the largest potential shift
        max_change = iterate();
        iteration++;
        
        // Print an update to the terminal every 100 sweeps so the user knows it hasn't frozen
        if (iteration % 100 == 0) {
            cout << "Iteration " << iteration << ", max change: " << max_change << endl;
        }
        
        // Failsafe break to prevent the program from running forever if it fails to converge
        if (iteration >= MAX_ITERATIONS) {
            cout << "Warning: Maximum iterations reached!\n";
            break;
        }
    // Keep looping as long as at least one node is changing by more than 0.000001 Volts
    } while (max_change > TOLERANCE);  
    
    cout << "Converged in " << iteration << " iterations.\n";
    return iteration; // Return the total sweep count for data logging
}

// Helper function to extract a 2D sheet of data from the 3D volume for Python/Gnuplot
void PotentialSolver::outputSlice(const string& filename, int plane_index, char axis) {
    ofstream file(filename); // Open a file stream to write to disk
    
    // Write metadata headers (lines starting with '#' are ignored by most plotting scripts)
    file << "# Electric potential slice\n";
    file << "# Axis: " << axis << " at index " << plane_index << endl;
    
    if (axis == 'z') {
        file << "# x(m)\ty(m)\tV(volts)\n";
        
        // Loop over the XY plane
        for (int i = 0; i < GRID_SIZE; i++) {
            for (int j = 0; j < GRID_SIZE; j++) {
                // Re-calculate the physical coordinates for the output text
                double x = GRID_MIN + i * GRID_SPACING;
                double y = GRID_MIN + j * GRID_SPACING;
                
                // Write a row: X_coord [tab] Y_coord [tab] Voltage at (x, y, fixed_z)
                file << x << "\t" << y << "\t" << V[i][j][plane_index] << "\n";
            }
            // Add a blank line after every row of X data (required formatting for plotting tools like pm3d)
            file << "\n";
        }
    }
    
    file.close(); // Save and close
    cout << "Output written to " << filename << endl;
}