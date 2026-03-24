/*
Stokes.cpp
-Solves the Navier-Stokes equations for CFD flow

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/


#include "Poisson.h"
#include "Stokes.h"
#include <iostream>

using namespace std;

Stokes::Stokes(const Mesh& mesh, double nu) : m_nu(nu) {
    m_nx = mesh.getNx();
    m_ny = mesh.getNy();
    int size = m_nx * m_ny;

    m_u.assign(size, 0.0);
    m_v.assign(size, 0.0);
    m_uStar.assign(size, 0.0);
    m_vStar.assign(size, 0.0);
}

void Stokes::predict(double dt, double dxi, double dyi) {
    double dxi2 = dxi * dxi;
    double dyi2 = dyi * dyi;

    // Loop through internal nodes to calculate uStar (u-momentum)
    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) {
            int k = i + j * m_nx;

            // 1. Convection terms (u * du/dx + v * du/dy)
            double convection = m_u[k] * (m_u[k+1] - m_u[k-1]) * 0.5 * dxi +
                                m_v[k] * (m_u[k+m_nx] - m_u[k-m_nx]) * 0.5 * dyi;

            // 2. Diffusion terms (nu * (d2u/dx2 + d2u/dy2))
            double diffusion = m_nu * (
                (m_u[k+1] - 2.0*m_u[k] + m_u[k-1]) * dxi2 +
                (m_u[k+m_nx] - 2.0*m_u[k] + m_u[k-m_nx]) * dyi2
            );

            m_uStar[k] = m_u[k] + dt * (-convection + diffusion);
        }
    }

    // Repeat similar logic for vStar (v-momentum)
    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) {
            int k = i + j * m_nx;

            double convection = m_u[k] * (m_v[k+1] - m_v[k-1]) * 0.5 * dxi +
                                m_v[k] * (m_v[k+m_nx] - m_v[k-m_nx]) * 0.5 * dyi;

            double diffusion = m_nu * (
                (m_v[k+1] - 2.0*m_v[k] + m_v[k-1]) * dxi2 +
                (m_v[k+m_nx] - 2.0*m_v[k] + m_v[k-m_nx]) * dyi2
            );

            m_vStar[k] = m_v[k] + dt * (-convection + diffusion);
        }
    }
}

void Stokes::correct(const std::vector<double>& p, double dt, double rho, double dxi, double dyi) {
    // Correct u-velocity
    for (int j = 0; j < m_ny; ++j) {
        for (int i = 1; i < m_nx; ++i) { // Skip i=0 for boundary
            int k = i + j * m_nx;
            // Gradient dp/dx uses pressure at k and k-1
            m_u[k] = m_uStar[k] - (dt / rho) * (p[k] - p[k-1]) * dxi;
        }
    }

    // Correct v-velocity
    for (int j = 1; j < m_ny; ++j) { // Skip j=0 for boundary
        for (int i = 0; i < m_nx; ++i) {
            int k = i + j * m_nx;
            // Gradient dp/dy uses pressure at k and k-nx
            m_v[k] = m_vStar[k] - (dt / rho) * (p[k] - p[k-m_nx]) * dyi;
        }
    }
}