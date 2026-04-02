/*
Poisson.cpp
-Solves the Poisson equations for CFD flow

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#include "Poisson.h"

// Constructor for the Poisson solver, takes a Mesh object to set up grid properties
Poisson::Poisson(const Mesh& mesh) {
    // Extract number of grid points in x and y directions from the mesh
    m_nx = mesh.getNx();
    m_ny = mesh.getNy();
    
    // Extract physical spacing between grid points in x and y
    m_dx = mesh.getDx();
    m_dy = mesh.getDy();

    // Initialize the pressure field vector with zeros, sized for the total number of cells
    m_p.assign(m_nx * m_ny, 0.0);
    
    // Initialize the Right-Hand Side (RHS) vector with zeros
    m_RHS.assign(m_nx * m_ny, 0.0);
    
    // Initialize the solid boundary mask vector with 'false' (assumes fluid by default)
    m_isSolid.assign(m_nx * m_ny, false);

    // Pre-calculate Cylinder Mask once
    // Define the center coordinates (cx, cy) and squared radius (R2) of the cylinder
    float cx = 0.3, cy = 0.49, R2 = 0.1 * 0.1;

    // Loop through every y-coordinate in the grid
    for (int j = 0; j < m_ny; ++j) {
        // Loop through every x-coordinate in the grid
        for (int i = 0; i < m_nx; ++i) {
            // Calculate the physical x and y positions of the current cell
            float x = i * m_dx;
            float y = j * m_dy;
            
            // Check if the current cell falls inside the circle using the circle equation: (x-h)^2 + (y-k)^2 <= R^2
            if (((x - cx) * (x - cx) + (y - cy) * (y - cy)) <= R2) {
                // If it is inside the circle, mark this cell as a solid
                m_isSolid[i + j * m_nx] = true;
            }
        }
    }
}

// Function to iteratively solve the pressure Poisson equation
void Poisson::solve(int max_Iterations, float tol, float dxi, float dyi) {
    // Pre-calculate squared inverse grid spacing for x and y
    float dxi2 = dxi * dxi;
    float dyi2 = dyi * dyi;
    
    // Pre-calculate the denominator term (beta) used in the central difference formula
    float beta = 2.0 * (dxi2 + dyi2);
    
    // Set the Successive Over-Relaxation (SOR) factor to accelerate convergence (1.0 < omega < 2.0)
    float omega = 1.8; 
    
    // Variable to track the maximum change in pressure during an iteration
    float localMaxError; 

    // Spawn a team of threads to run the following block in parallel
    #pragma omp parallel
    {
        // Loop through the solver iterations up to the maximum allowed
        for (int iter = 0; iter < max_Iterations; ++iter) {
            
            // Ensure only one thread resets the error variable at the start of the iteration
            #pragma omp single
            {
                localMaxError = 0.0;
            }

            // PASS 1: RED CELLS (i+j is even)
            // Distribute the loop among threads and find the maximum error across all threads
            #pragma omp for reduction(max:localMaxError)
            for (int j = 1; j < m_ny - 1; ++j) {
                for (int i = 1; i < m_nx - 1; ++i) {
                    // Check if the cell belongs to the "Red" set (checkerboard pattern)
                    if ((i + j) % 2 == 0) {
                        // Calculate the 1D index for the current 2D (i, j) coordinates
                        int k = i + j * m_nx;
                        
                        // If the cell is inside the cylinder mask, force pressure to 0
                        if (m_isSolid[k]) {
                            m_p[k] = 0.0; 
                        } else {
                            // Store the old pressure value to calculate the error later
                            float p_old = m_p[k];
                            
                            // Sum the pressure of the 4 adjacent neighbors multiplied by their respective inverse spacing squared
                            float neighbors = dxi2 * (m_p[k+1] + m_p[k-1]) + 
                                               dyi2 * (m_p[k+m_nx] + m_p[k-m_nx]);
                                               
                            // Calculate the new Gauss-Seidel pressure estimate
                            float p_gs = (neighbors - m_RHS[k]) / beta;
                            
                            // Apply the SOR formula to blend the old pressure and the new estimate
                            m_p[k] = (1.0 - omega) * p_old + omega * p_gs;
                            
                            // Update the local max error if the change in this cell is the largest seen so far
                            localMaxError = std::max(localMaxError, std::abs(m_p[k] - p_old));
                        }
                    }
                }
            }

            // PASS 2: BLACK CELLS (i+j is odd)
            // Distribute the loop among threads and continue finding the maximum error
            #pragma omp for reduction(max:localMaxError)
            for (int j = 1; j < m_ny - 1; ++j) {
                for (int i = 1; i < m_nx - 1; ++i) {
                    // Check if the cell belongs to the "Black" set
                    if ((i + j) % 2 != 0) {
                        // Calculate the 1D index
                        int k = i + j * m_nx;
                        
                        // If the cell is solid, force pressure to 0
                        if (m_isSolid[k]) {
                            m_p[k] = 0.0;
                        } else {
                            // Store the old pressure value
                            float p_old = m_p[k];
                            
                            // Sum the neighbors using the newly updated "Red" values
                            float neighbors = dxi2 * (m_p[k+1] + m_p[k-1]) + 
                                               dyi2 * (m_p[k+m_nx] + m_p[k-m_nx]);
                                               
                            // Calculate the Gauss-Seidel estimate
                            float p_gs = (neighbors - m_RHS[k]) / beta;
                            
                            // Apply the SOR blend
                            m_p[k] = (1.0 - omega) * p_old + omega * p_gs;
                            
                            // Update the maximum error
                            localMaxError = std::max(localMaxError, std::abs(m_p[k] - p_old));
                        }
                    }
                }
            }

            // Apply Boundary Conditions & Check Convergence
            // Declare a static boolean to signal threads to break (persists across loop iterations)
            static bool shouldBreak;
            
            // Ensure only one thread executes the boundary condition assignments
            #pragma omp single
            {
                // Neumann BCs: dp/dn = 0 (Pressure at wall = Pressure just inside wall)
                // Loop through the x-axis to set Top and Bottom boundaries
                for (int i = 0; i < m_nx; ++i) {
                    m_p[i + 0 * m_nx] = m_p[i + 1 * m_nx];               // Bottom boundary copies the row above it
                    m_p[i + (m_ny-1) * m_nx] = m_p[i + (m_ny-2) * m_nx]; // Top boundary copies the row below it
                }
                
                // Loop through the y-axis to set Inlet and Outlet boundaries
                for (int j = 0; j < m_ny; ++j) {
                    m_p[0 + j * m_nx] = m_p[1 + j * m_nx];               // Inlet boundary copies the column to its right
                    m_p[(m_nx-1) + j * m_nx] = 0.0;                      // Outlet boundary is fixed to 0 (Dirichlet BC)
                }
            }
            
            // Check if the maximum error across the grid is less than the specified tolerance
            // If it is, break out of the iteration loop early
            if (localMaxError < tol) break;
            
            // Force all threads to wait here until everyone has reached this point
            #pragma omp barrier 
            
            // If the static break flag was set (though not explicitly set to true in this snippet), break the loop
            if (shouldBreak) break;
        }
    }
}

// Function to build the Right-Hand Side of the Poisson equation based on intermediate velocities
void Poisson::buildRHS(const std::vector<float>& uStar, 
                       const std::vector<float>& vStar, 
                       float dt, float rho, float dxi, float dyi) 
{
    // Reset the entire RHS vector to zero before calculating new values
    std::fill(m_RHS.begin(), m_RHS.end(), 0.0);

    // Loop through the inner grid cells (ignoring the outer boundary layer)
    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) {
            // Calculate the 1D index for the current cell
            int k = i + j * m_nx;

            // Correct 1D indices for the neighboring cells
            int k_right = (i + 1) + j * m_nx; // Right neighbor
            int k_left  = (i - 1) + j * m_nx; // Left neighbor
            int k_up    = i + (j + 1) * m_nx; // Top neighbor
            int k_down  = i + (j - 1) * m_nx; // Bottom neighbor

            // Calculate the horizontal velocity gradient (du/dx) using central differencing
            float du_dx = (uStar[k_right] - uStar[k_left]) * 0.5 * dxi;
            
            // Calculate the vertical velocity gradient (dv/dy) using central differencing
            float dv_dy = (vStar[k_up] - vStar[k_down]) * 0.5 * dyi;

            // Calculate the source term for the pressure Poisson equation
            // Notice: Positive sign here balances the subtraction in the solver!
            // It represents the divergence of the intermediate velocity field scaled by density and time step
            m_RHS[k] = (rho / dt) * (du_dx + dv_dy);
        }
    }
}