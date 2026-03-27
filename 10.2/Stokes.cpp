/*
Stokes.cpp
-Solves the Navier-Stokes equations for CFD flow

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/


#include "Poisson.h"
#include "Stokes.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

Stokes::Stokes(int nx, int ny, int imin, int jmin, 
               int imax, int jmax, const Mesh& mesh, double nu) 
    : m_nx(nx), m_ny(ny)
      m_imin(imin), m_jmin(jmin),
      m_imax(imax), m_jmax(jmax),
      m_nu(nu)
{

}

void Stokes::predict(double dt, double dxi, double dyi) {


    int size = m_nx * m_ny;

    // U-grid
    m_u.assign(size, 0.0);
    m_uStar.assign(size, 0.0);

    // V-grid
    m_v.assign(size, 0.0);
    m_vStar.assign(size, 0.0);

    double dxi2 = dxi * dxi;
    double dyi2 = dyi * dyi;

    // U-momentum: Solve for uStar
    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) { 
            int k = i + j * m_nx;

            double v_here = 0.25 * (m_v[k-1] + m_v[k-1+m_nx] + m_v[k] + m_v[k+m_nx]);

            double u_grad_x = (m_u[k] > 0) ? (m_u[k] - m_u[k-1]) * dxi : (m_u[k+1] - m_u[k]) * dxi;
            double u_grad_y = (v_here > 0) ? (m_u[k] - m_u[k-m_nx]) * dyi : (m_u[k+m_nx] - m_u[k]) * dyi;
            double convection = m_u[k] * u_grad_x + v_here * u_grad_y;
            
            double diffusion = m_nu * (
                (m_u[k+1] - 2.0*m_u[k] + m_u[k-1]) * dxi2 +
                (m_u[k+m_nx] - 2.0*m_u[k] + m_u[k-m_nx]) * dyi2
            );

            m_uStar[k] = m_u[k] + dt * (-convection + diffusion);
        }
    }

    // V-momentum: Solve for vStar
    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) {
            int k = i + j * m_nx;

            double u_here = 0.25 * (m_u[k-m_nx] + m_u[k] + m_u[k+1-m_nx] + m_u[k+1]);
            double v_grad_x = (m_v[k] > 0) ? (m_v[k] - m_v[k-1]) * dxi : (m_v[k+1] - m_v[k]) * dxi;
            double v_grad_y = (u_here > 0) ? (m_v[k] - m_v[k-m_nx]) * dyi : (m_v[k+m_nx] - m_v[k]) * dyi;
            double convection = m_v[k] * v_grad_x + u_here * v_grad_y;
            double diffusion = m_nu * (
                (m_v[k+1] - 2.0*m_v[k] + m_v[k-1]) * dxi2 +
                (m_v[k+m_nx] - 2.0*m_v[k] + m_v[k-m_nx]) * dyi2
            );

            m_vStar[k] = m_v[k] + dt * (-convection + diffusion);
        }
    }
}

void Stokes::applyBoundary(double U_lid){
// 1. Bottom Wall (j = 0)
    for (int i = 0; i < m_nx; ++i) {
        int k = i + 0 * m_nx;
        m_u[k] = 0.0;
        m_v[k] = 0.0;
    }

    // 2. Top Wall (j = ny - 1)
    for (int i = 0; i < m_nx; ++i) {
        int k = i + (m_ny - 1) * m_nx;
        m_u[k] = U_lid;
        m_v[k] = 0.0;
    }

    // 3. Left Wall (i = 0)
    for (int j = 0; j < m_ny; ++j) {
        int k = 0 + j * m_nx;
        m_u[k] = 0.0;
        m_v[k] = 0.0;
    }

    // 4. Right Wall (i = nx - 1)
    for (int j = 0; j < m_ny; ++j) {
        int k = (m_nx - 1) + j * m_nx;
        m_u[k] = 0.0;
        m_v[k] = 0.0;
    }
}

void Stokes::applyBoundaryToStar(double U_lid) {
    // Bottom (j=0): no-slip
    for (int i = 0; i < m_nx; ++i) {
        m_uStar[i] = 0.0;
        m_vStar[i] = 0.0;
    }
    // Top (j=ny-1): moving lid
    for (int i = 0; i < m_nx; ++i) {
        int k = i + (m_ny - 1) * m_nx;
        m_uStar[k] = U_lid;
        m_vStar[k] = 0.0;
    }
    // Left (i=0): no-slip
    for (int j = 0; j < m_ny; ++j) {
        m_uStar[j * m_nx] = 0.0;
        m_vStar[j * m_nx] = 0.0;
    }
    // Right (i=nx-1): no-slip
    for (int j = 0; j < m_ny; ++j) {
        int k = (m_nx - 1) + j * m_nx;
        m_uStar[k] = 0.0;
        m_vStar[k] = 0.0;
    }
}

void Stokes::correct(const std::vector<double>& p, double dt, double rho, double dxi, double dyi) {
    int n = 1;
    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 2; i < m_nx - 1; ++i) {
            n += 1;
            m_u[n] = m_uStar[i,j] - (dt / rho) * (p[i,j] - p[i-1, j]) * dxi;
        }
    }
    int n = 1;
    for (int j = 2; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) {
            n += 1;
            m_u[n] = m_uStar[i,j] - (dt / rho) * (p[i,j] - p[i, j-1]) * dxi;
        }
    }
}

void Stokes::exportFrame(std::string filename, int frame, const std::vector<double>& p) {
    std::ofstream file(filename);
    file << "x,y,u,v,p\n";

    double dx = 1.0 / m_nx;
    double dy = 1.0 / m_ny;

    for (int j = 0; j < m_ny; ++j) {
        for (int i = 0; i < m_nx; ++i) {
            int k = i + j * m_nx;
            file << i * dx << "," << j * dy << "," 
                 << m_u[k] << "," << m_v[k] << "," << p[k] << "\n";
        }
    }
    file.close();
}
