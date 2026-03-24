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
    m_nx = mesh.getNx();
    m_ny = mesh.getNy();
    m_size = m_nx*m_ny;

    m_p.assign(m_size, 0.0);
    m_RHS.assign(m_size, 0.0);
    m_diag.resize(m_size);

    double dxi2 = mesh.getDxi() * mesh.getDxi();
    double dyi2 = mesh.getDyi() * mesh.getDyi();

    m_rowPtr.push_back(0); // First row starts at index 0

    for (int j = 0; j < m_ny; ++j) {
        for (int i = 0; i < m_nx; ++i) {
            int k = i + j * m_nx;

            // Calculating the center coefficents with Neumann Adjusments
            double center = 2.0 * (dxi2 + dyi2);
            if (j == 0) center -= dyi2;        // Bottom Wall
            if (j == m_ny - 1) center -= dyi2; // Top Wall
            if (i == 0) center -= dxi2;        // Left Wall
            if (i == m_nx - 1) center -= dxi2; // Right Wall
            
            // Adding neighbors in increasing column order (this is important)

            // 1. Bottom neighbor (j-1)
            if (j > 0) {
                m_values.push_back(dyi2);
                m_colIndices.push_back(k - m_nx);
            }
            
            // 2. Left neighbor (i-1)
            if (i > 0) {
                m_values.push_back(dxi2);
                m_colIndices.push_back(k - 1);
            }

            // 3. Center Node
            m_values.push_back(center);
            m_colIndices.push_back(k);
            m_diag[k] = center;

            // 4. Right neighbor (i+1)
            if (i < m_nx - 1) {
                m_values.push_back(dxi2);
                m_colIndices.push_back(k + 1);
            }

            // 5. Top neighbor (j+1)
            if (j < m_ny - 1) {
                m_values.push_back(dyi2);
                m_colIndices.push_back(k + m_nx);
            }

            // Mark where the next row starts
            m_rowPtr.push_back(m_values.size());
        }
    }
    m_values[m_rowPtr[0]] = 1.0;
}

void Poisson::buildRHS(const std::vector<double>& uStar, 
                       const std::vector<double>& vStar, 
                       double dt, double rho, double dxi, double dyi) 
{
    for (int j = 0; j < m_ny; ++j) {
        for (int i = 0; i < m_nx; ++i) {
            int k = i + j * m_nx;

            // In a staggered grid:
            // uStar[k] is the velocity on the RIGHT face of cell k
            // uStar[k-1] is the velocity on the LEFT face of cell k
            double du_dx = (i > 0) ? (uStar[k] - uStar[k-1]) * dxi : uStar[k] * dxi;
            double dv_dy = (j > 0) ? (vStar[k] - vStar[k-m_nx]) * dyi : vStar[k] * dyi;

            m_RHS[k] = (rho / dt) * (du_dx + dv_dy);
        }
    }
    m_RHS[0] = 0.0;
}