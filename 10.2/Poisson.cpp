/*
Poisson.cpp
-Solves the Poisson equations for CFD flow

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/


#include "Poisson.h"
#include <iostream>

using namespace std;

Poisson::Poisson(const Mesh& mesh){

    m_dxi2 = mesh.getDxi() * mesh.getDxi();
    m_dyi2 = mesh.getDyi() * mesh.getDyi();
    m_nx = mesh.getNx();
    m_ny = mesh.getNy();
    m_size = m_nx*m_ny;

    m_p.assign(m_size, 0.0);
    m_RHS.assign(m_size, 0.0);
    m_diag.resize(m_size);

    m_rowPtr.push_back(0); // First row starts at index 0

    for (int j = 0; j < m_ny; ++j) {
        for (int i = 0; i < m_nx; ++i) {
            int k = i + j * m_nx;

            double dxi2 = m_dxi2;
            double dyi2 = m_dyi2;
            // Calculating the center coefficents with Neumann Adjusments
            double center = 2.0 * (dxi2 + dyi2);
            if (j == 0) center -= dyi2;        // Bottom Wall
            if (j == m_ny - 1) center -= dyi2; // Top Wall
            if (i == 0) center -= dxi2;        // Left Wall
            if (i == m_nx - 1) center -= dxi2; // Right Wall
            
            // Adding neighbors in increasing column order (this is important)

            // 1. Bottom neighbor (j-1)
            if (j > 0) {
                m_values.push_back(-dyi2);
                m_colIndices.push_back(k - m_nx);
            }
            
            // 2. Left neighbor (i-1)
            if (i > 0) {
                m_values.push_back(-dxi2);
                m_colIndices.push_back(k - 1);
            }

            // 3. Center Node
            m_values.push_back(center);
            m_colIndices.push_back(k);
            m_diag[k] = center;

            // 4. Right neighbor (i+1)
            if (i < m_nx - 1) {
                m_values.push_back(-dxi2);
                m_colIndices.push_back(k + 1);
            }

            // 5. Top neighbor (j+1)
            if (j < m_ny - 1) {
                m_values.push_back(-dyi2);
                m_colIndices.push_back(k + m_nx);
            }

            // Mark where the next row starts
            m_rowPtr.push_back(m_values.size());
        }
    }
    m_values[m_rowPtr[0]] = 1.0;
}

void Poisson::solve(int max_Iterations, double tol){

    double dxi2 = m_dxi2;
    double dyi2 = m_dyi2;
    double beta = 2.0 * (dxi2 + dyi2);

    for (int iter = 0; iter < max_Iterations; ++iter) {
        double maxError = 0.0;
    

        // Loop through internal cells (Subtraction Method safe zone)
        for (int j = 1; j < m_ny - 1; ++j) {
            for (int i = 1; i < m_nx - 1; ++i) {
                int k = i + j * m_nx;
                
                double oldP = m_p[k];
                
                // Gauss-Seidel Formula
                double neighbors = dxi2 * (m_p[k+1] + m_p[k-1]) + 
                                   dyi2 * (m_p[k+m_nx] + m_p[k-m_nx]);
                
                m_p[k] = (neighbors - m_RHS[k]) / beta;

                maxError = std::max(maxError, std::abs(m_p[k] - oldP));
            }
        }
        // Including a copy for the boundaries
        for (int i = 0; i < m_nx; ++i) {
            m_p[i + 0 * m_nx] = m_p[i + 1 * m_nx];                   // Bottom
            m_p[i + (m_ny-1) * m_nx] = m_p[i + (m_ny-2) * m_nx];     // Top
        }
        for (int j = 0; j < m_ny; ++j) {
            m_p[0 + j * m_nx] = m_p[1 + j * m_nx];                   // Left
            m_p[(m_nx-1) + j * m_nx] = m_p[(m_nx-2) + j * m_nx];     // Right
        }

        //m_p[0] = 0.0;

        // Check for convergence
        if (maxError < tol) break;
    }
}

void Poisson::buildRHS(const std::vector<double>& uStar, 
                       const std::vector<double>& vStar, 
                       double dt, double rho, double dxi, double dyi) 
{
    fill(m_RHS.begin(), m_RHS.end(), 0.0);
    int n = 0;
    for (int j = 1; j < m_ny-1; ++j) {
        for (int i = 1; i < m_nx-1; ++i) {
            //int k = i + j * m_nx;
            n += 1;
            double du_dx = (uStar[i+1, j] - uStar[i,j]) * dxi;
            double dv_dy = (vStar[i, j+1] - vStar[i, j]) * dyi;

            m_RHS[n] = -(rho / dt) * (du_dx + dv_dy);


        }
    }
}