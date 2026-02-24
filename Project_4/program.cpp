#include "program.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

// Constructor: Initialize the potential and fixed arrays, and define point charges
PotentialSolver::PotentialSolver() {
    // Initialize 3D potential array with all zeros
    V.resize(GRID_SIZE, vector<vector<double>>(GRID_SIZE, vector<double>(GRID_SIZE, 0.0)));
    
    // Initialize 3D fixed array with all false (no points fixed initially)
    fixed.resize(GRID_SIZE, vector<vector<bool>>(GRID_SIZE, vector<bool>(GRID_SIZE, false)));
    
    // Define the four point charges in the system
    charges.push_back({0.25, 0.0, 0.0, CHARGE});    // +1 uC at (25,0,0) cm = (0.25,0,0) m
    charges.push_back({0.0, 0.25, 0.0, CHARGE});    // +1 uC at (0,25,0) cm = (0,0.25,0) m
    charges.push_back({-0.25, 0.0, 0.0, -CHARGE});  // -1 uC at (-25,0,0) cm = (-0.25,0,0) m
    charges.push_back({0.0, -0.25, 0.0, -CHARGE});  // -1 uC at (0,-25,0) cm = (0,-0.25,0) m
}

// Set up boundary conditions and initial potential distribution
void PotentialSolver::initializeConditions() {
    // Set boundary conditions: all six faces of the cube held at zero potential
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            // Fix potential at x = GRID_MIN face
            fixed[0][i][j] = true;
            // Fix potential at x = GRID_MAX face
            fixed[GRID_SIZE-1][i][j] = true;
            // Fix potential at y = GRID_MIN face
            fixed[i][0][j] = true;
            // Fix potential at y = GRID_MAX face
            fixed[i][GRID_SIZE-1][j] = true;
            // Fix potential at z = GRID_MIN face
            fixed[i][j][0] = true;
            // Fix potential at z = GRID_MAX face
            fixed[i][j][GRID_SIZE-1] = true;
        }
    }
    
    // Set up the conducting disk at z=0 plane with radius 10 cm at 0.25 V
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            // Convert grid indices to physical coordinates
            double x = GRID_MIN + i * GRID_SPACING;
            double y = GRID_MIN + j * GRID_SPACING;
            
            // Calculate distance from z-axis (radial distance in xy-plane)
            double r = sqrt(x*x + y*y);
            
            // Index for z = 0 plane (middle of the grid)
            int k = GRID_SIZE / 2;
            
            // If point is within disk radius, set voltage and mark as fixed
            if (r <= DISK_RADIUS) {
                V[i][j][k] = DISK_VOLTAGE;  // Set disk voltage
                fixed[i][j][k] = true;       // Mark as fixed boundary condition
            }
        }
    }
    
    // Initialize potential due to point charges throughout the grid
    for (int i = 1; i < GRID_SIZE-1; i++) {           // Skip boundary points
        for (int j = 1; j < GRID_SIZE-1; j++) {       // Skip boundary points
            for (int k = 1; k < GRID_SIZE-1; k++) {   // Skip boundary points
                // Only set initial potential at non-fixed points
                if (!fixed[i][j][k]) {
                    // Convert grid indices to physical coordinates
                    double x = GRID_MIN + i * GRID_SPACING;
                    double y = GRID_MIN + j * GRID_SPACING;
                    double z = GRID_MIN + k * GRID_SPACING;
                    
                    // Add contribution from each point charge using V = kq/r
                    for (const auto& charge : charges) {
                        // Calculate displacement from charge to grid point
                        double dx = x - charge.x;
                        double dy = y - charge.y;
                        double dz = z - charge.z;
                        
                        // Calculate distance from charge to grid point
                        double r = sqrt(dx*dx + dy*dy + dz*dz);
                        
                        // Avoid singularity at charge location
                        if (r > GRID_SPACING/2) {
                            // Add potential contribution: V = kq/r
                            V[i][j][k] += K * charge.q / r;
                        }
                    }
                }
            }
        }
    }
}

// Perform one iteration of the successive over-relaxation algorithm
double PotentialSolver::iterate() {
    double max_change = 0.0;  // Track maximum potential change in this iteration
    
    // Loop through all interior grid points (skip boundaries)
    for (int i = 1; i < GRID_SIZE-1; i++) {
        for (int j = 1; j < GRID_SIZE-1; j++) {
            for (int k = 1; k < GRID_SIZE-1; k++) {
                // Only update non-fixed points
                if (!fixed[i][j][k]) {
                    // Store old potential value
                    double old_V = V[i][j][k];
                    
                    // Calculate new potential as average of six neighbors
                    // This satisfies Laplace's equation: ∇²V = 0
                    double V_new = (V[i+1][j][k] + V[i-1][j][k] +  // x-neighbors
                                   V[i][j+1][k] + V[i][j-1][k] +    // y-neighbors
                                   V[i][j][k+1] + V[i][j][k-1]) / 6.0; // z-neighbors
                    
                    // Apply over-relaxation: V = V_old + ω(V_new - V_old)
                    // This accelerates convergence compared to simple relaxation
                    V[i][j][k] = old_V + OMEGA * (V_new - old_V);
                    
                    // Calculate absolute change in potential
                    double change = fabs(V[i][j][k] - old_V);
                    
                    // Update maximum change if this point changed more
                    if (change > max_change) {
                        max_change = change;
                    }
                }
            }
        }
    }
    
    // Return the maximum change for convergence checking
    return max_change;
}

// Solve for potential distribution until convergence
void PotentialSolver::solve() {
    // Print solver configuration
    cout << "Starting over-relaxation solver...\n";
    cout << "Grid size: " << GRID_SIZE << "x" << GRID_SIZE << "x" << GRID_SIZE << endl;
    cout << "Omega: " << OMEGA << endl;
    
    int iteration = 0;      // Iteration counter
    double max_change;      // Maximum potential change in current iteration
    
    // Iterate until convergence or maximum iterations reached
    do {
        // Perform one iteration and get maximum change
        max_change = iterate();
        iteration++;
        
        // Print progress every 100 iterations
        if (iteration % 100 == 0) {
            cout << "Iteration " << iteration << ", max change: " << max_change << endl;
        }
        
        // Check if maximum iterations exceeded
        if (iteration >= MAX_ITERATIONS) {
            cout << "Warning: Maximum iterations reached!\n";
            break;
        }
    } while (max_change > TOLERANCE);  // Continue until changes are small enough
    
    cout << "Converged in " << iteration << " iterations.\n";
}

// Output a 2D slice of the potential to a data file
void PotentialSolver::outputSlice(const string& filename, int plane_index, char axis) {
    ofstream file(filename);  // Open output file
    
    // Write header information
    file << "# Electric potential slice\n";
    file << "# Axis: " << axis << " at index " << plane_index << endl;
    
    // Currently only supports z-axis slices
    if (axis == 'z') {
        file << "# x(m)\ty(m)\tV(volts)\n";
        
        // Loop through x and y coordinates
        for (int i = 0; i < GRID_SIZE; i++) {
            for (int j = 0; j < GRID_SIZE; j++) {
                // Convert indices to physical coordinates
                double x = GRID_MIN + i * GRID_SPACING;
                double y = GRID_MIN + j * GRID_SPACING;
                
                // Write x, y, and potential value
                file << x << "\t" << y << "\t" << V[i][j][plane_index] << "\n";
            }
            // Blank line between rows for gnuplot pm3d format
            file << "\n";
        }
    }
    
    file.close();
    cout << "Output written to " << filename << endl;
}

// Output full 3D potential data in VTK format for ParaView
void PotentialSolver::outputVTK(const string& filename) {
    ofstream file(filename);  // Open output file
    
    // Write VTK file header
    file << "# vtk DataFile Version 3.0\n";
    file << "Electric Potential\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    
    // Define grid dimensions
    file << "DIMENSIONS " << GRID_SIZE << " " << GRID_SIZE << " " << GRID_SIZE << "\n";
    
    // Define grid origin (minimum corner)
    file << "ORIGIN " << GRID_MIN << " " << GRID_MIN << " " << GRID_MIN << "\n";
    
    // Define grid spacing
    file << "SPACING " << GRID_SPACING << " " << GRID_SPACING << " " << GRID_SPACING << "\n";
    
    // Define scalar field data
    file << "POINT_DATA " << (GRID_SIZE * GRID_SIZE * GRID_SIZE) << "\n";
    file << "SCALARS potential double 1\n";
    file << "LOOKUP_TABLE default\n";
    
    // Write potential values in k-j-i order (z-y-x) as required by VTK
    for (int k = 0; k < GRID_SIZE; k++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            for (int i = 0; i < GRID_SIZE; i++) {
                file << V[i][j][k] << "\n";
            }
        }
    }
    
    file.close();
    cout << "VTK output written to " << filename << endl;
}

// Output full 3D potential data in gnuplot-compatible format
void PotentialSolver::outputGnuplot3D(const string& filename) {
    ofstream file(filename);  // Open output file

    
    
    // Write header
    file << "# 3D Electric potential data for gnuplot\n";
    file << "# x(m)\ty(m)\tz(m)\tV(volts)\n";
    
    // Loop through entire 3D grid
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            for (int k = 0; k < GRID_SIZE; k++) {
                // Convert indices to physical coordinates
                double x = GRID_MIN + i * GRID_SPACING;
                double y = GRID_MIN + j * GRID_SPACING;
                double z = GRID_MIN + k * GRID_SPACING;
                
                // Write x, y, z, and potential value
                file << x << "\t" << y << "\t" << z << "\t" << V[i][j][k] << "\n";
            }
        }
    }
    
    file.close();
    cout << "Gnuplot 3D output written to " << filename << endl;
}