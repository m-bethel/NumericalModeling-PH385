/*
Poisson.cpp
-Solves the Poisson equations for CFD flow

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/


#include "Poisson.h"

Poisson::Poisson(const Mesh& mesh) {
    m_nx = mesh.getNx();
    m_ny = mesh.getNy();
    m_dx = mesh.getDx();
    m_dy = mesh.getDy();

    m_p.assign(m_nx * m_ny, 0.0);
    m_RHS.assign(m_nx * m_ny, 0.0);
    m_isSolid.assign(m_nx * m_ny, false);

    // Pre-calculate Cylinder Mask once
    float cx = 0.3, cy = 0.49, R2 = 0.1 * 0.1;

    for (int j = 0; j < m_ny; ++j) {
        for (int i = 0; i < m_nx; ++i) {
            float x = i * m_dx;
            float y = j * m_dy;
            if (((x - cx) * (x - cx) + (y - cy) * (y - cy)) <= R2) {
                m_isSolid[i + j * m_nx] = true;
            }
        }
    }
}

void Poisson::solve(int max_Iterations, float tol, float dxi, float dyi) {
    float dxi2 = dxi * dxi;
    float dyi2 = dyi * dyi;
    float beta = 2.0 * (dxi2 + dyi2);
    float omega = 1.8; // Successive Over-Relaxation factor
    float localMaxError; 

    #pragma omp parallel
    {
        for (int iter = 0; iter < max_Iterations; ++iter) {
            #pragma omp single
            {
                localMaxError = 0.0;
            }

            // PASS 1: RED CELLS (i+j is even)
            #pragma omp for reduction(max:localMaxError)
            for (int j = 1; j < m_ny - 1; ++j) {
                for (int i = 1; i < m_nx - 1; ++i) {
                    if ((i + j) % 2 == 0) {
                        int k = i + j * m_nx;
                        if (m_isSolid[k]) {
                            m_p[k] = 0.0; 
                        } else {
                            float p_old = m_p[k];
                            float neighbors = dxi2 * (m_p[k+1] + m_p[k-1]) + 
                                               dyi2 * (m_p[k+m_nx] + m_p[k-m_nx]);
                            float p_gs = (neighbors - m_RHS[k]) / beta;
                            m_p[k] = (1.0 - omega) * p_old + omega * p_gs;
                            localMaxError = std::max(localMaxError, std::abs(m_p[k] - p_old));
                        }
                    }
                }
            }

            // PASS 2: BLACK CELLS (i+j is odd)
            #pragma omp for reduction(max:localMaxError)
            for (int j = 1; j < m_ny - 1; ++j) {
                for (int i = 1; i < m_nx - 1; ++i) {
                    if ((i + j) % 2 != 0) {
                        int k = i + j * m_nx;
                        if (m_isSolid[k]) {
                            m_p[k] = 0.0;
                        } else {
                            float p_old = m_p[k];
                            float neighbors = dxi2 * (m_p[k+1] + m_p[k-1]) + 
                                               dyi2 * (m_p[k+m_nx] + m_p[k-m_nx]);
                            float p_gs = (neighbors - m_RHS[k]) / beta;
                            m_p[k] = (1.0 - omega) * p_old + omega * p_gs;
                            localMaxError = std::max(localMaxError, std::abs(m_p[k] - p_old));
                        }
                    }
                }
            }

            // Apply Boundary Conditions & Check Convergence
            static bool shouldBreak;
            #pragma omp single
            {
                // Neumann BCs: dp/dn = 0 (Pressure at wall = Pressure just inside wall)
                for (int i = 0; i < m_nx; ++i) {
                    m_p[i + 0 * m_nx] = m_p[i + 1 * m_nx];               // Bottom
                    m_p[i + (m_ny-1) * m_nx] = m_p[i + (m_ny-2) * m_nx]; // Top
                }
                for (int j = 0; j < m_ny; ++j) {
                    m_p[0 + j * m_nx] = m_p[1 + j * m_nx];               // Inlet
                    m_p[(m_nx-1) + j * m_nx] = 0.0;                      // Outlet Fixed
                }
            }
            
            // Check if we should stop early
            if (localMaxError < tol) break;
            #pragma omp barrier 
            if (shouldBreak) break;
        }
    }
}

void Poisson::buildRHS(const std::vector<float>& uStar, 
                       const std::vector<float>& vStar, 
                       float dt, float rho, float dxi, float dyi) 
{
    // Reset RHS to zero
    std::fill(m_RHS.begin(), m_RHS.end(), 0.0);

    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) {
            int k = i + j * m_nx;

            // Correct 1D indices for neighbors
            int k_right = (i + 1) + j * m_nx;
            int k_left  = (i - 1) + j * m_nx;
            int k_up    = i + (j + 1) * m_nx;
            int k_down  = i + (j - 1) * m_nx;

            // Central difference calculation for divergence
            float du_dx = (uStar[k_right] - uStar[k_left]) * 0.5 * dxi;
            float dv_dy = (vStar[k_up] - vStar[k_down]) * 0.5 * dyi;

            // Notice: Positve sign here balances the subtraction in the solver!
            m_RHS[k] = (rho / dt) * (du_dx + dv_dy);
        }
    }
}