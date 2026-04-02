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

Stokes::Stokes(const Mesh& mesh, float nu) : m_nu(nu) {
    m_nx = mesh.getNx();
    m_ny = mesh.getNy();

    int size = m_nx * m_ny;

    // U-grid
    m_u.assign(size, 0.0);
    m_uStar.assign(size, 0.0);

    // V-grid
    m_v.assign(size, 0.0);
    m_vStar.assign(size, 0.0);
}

void Stokes::predict(float dt, float dxi, float dyi) {
    float dxi2 = dxi * dxi;
    float dyi2 = dyi * dyi;
    float dx = 1.0 / dxi;
    float dy = 1.0 / dyi;

    // Cylinder Parameters (Match these in your correct() function)
    float cx = 0.3;  // Center X
    float cy = 0.49f;  // Center Y
    float R = 0.1;   // Radius
    float R2 = R * R;

    // U-momentum: Solve for uStar
    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) { 
            int k = i + j * m_nx;

            // Physical coordinates for the mask
            float x = i * dx; 
            float y = (j + 0.5) * dy;

            if (((x - cx) * (x - cx) + (y - cy) * (y - cy)) <= R2) {
                m_uStar[k] = 0.0;
            } else {
                // Standard Advection + Diffusion
                float v_here = 0.25 * (m_v[k] + m_v[k + 1] + m_v[k - m_nx] + m_v[k + 1 - m_nx]);
                float u_grad_x = (m_u[k] > 0) ? (m_u[k] - m_u[k-1]) * dxi : (m_u[k+1] - m_u[k]) * dxi;
                float u_grad_y = (v_here > 0) ? (m_u[k] - m_u[k-m_nx]) * dyi : (m_u[k+m_nx] - m_u[k]) * dyi;
                
                float convection = m_u[k] * u_grad_x + v_here * u_grad_y;
                float diffusion = m_nu * (
                    (m_u[k+1] - 2.0*m_u[k] + m_u[k-1]) * dxi2 +
                    (m_u[k+m_nx] - 2.0*m_u[k] + m_u[k-m_nx]) * dyi2
                );

                m_uStar[k] = m_u[k] + dt * (-convection + diffusion);
            }
        }
    }

    // V-momentum: Solve for vStar
    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) {
            int k = i + j * m_nx;

            // Physical coordinates for the mask
            float x = (i + 0.5) * dx;
            float y = j * dy;

            if (((x - cx) * (x - cx) + (y - cy) * (y - cy)) <= R2) {
                m_vStar[k] = 0.0;
            } else {
                float u_here = 0.25 * (m_u[k] + m_u[k - 1] + m_u[k + m_nx] + m_u[k + m_nx - 1]);
                float v_here = m_v[k];
                float v_grad_x = (u_here > 0) ? (m_v[k] - m_v[k-1]) * dxi : (m_v[k+1] - m_v[k]) * dxi;
                float v_grad_y = (m_v[k] > 0) ? (m_v[k] - m_v[k-m_nx]) * dyi : (m_v[k+m_nx] - m_v[k]) * dyi;
                
                float convection = v_here * v_grad_y + u_here * v_grad_x; // Corrected order
                float diffusion = m_nu * (
                    (m_v[k+1] - 2.0*m_v[k] + m_v[k-1]) * dxi2 +
                    (m_v[k+m_nx] - 2.0*m_v[k] + m_v[k-m_nx]) * dyi2
                );

                m_vStar[k] = m_v[k] + dt * (-convection + diffusion);
            }
        }
    }
}

void Stokes::applyBoundary(float U_in){
// Bottom (j=0) and Top (j=ny-1): Solid walls (No-slip)
    for (int i = 0; i < m_nx; ++i) {
        m_u[i + 0 * m_nx] = 0.0;                //  m_u instead of m_uStar
        m_v[i + 0 * m_nx] = 0.0;           
        m_u[i + (m_ny - 1) * m_nx] = 0.0;       
        m_v[i + (m_ny - 1) * m_nx] = 0.0;
    }

    // Left (i=0): Inlet (Fluid entering)
    for (int j = 0; j < m_ny; ++j) {
        m_u[0 + j * m_nx] = U_in;               
        m_v[0 + j * m_nx] = 0.0;
    }
    
    // Right (i=nx-1): Outlet (Zero-gradient / free exit)
    for (int j = 0; j < m_ny; ++j) {
        int k_out = (m_nx - 1) + j * m_nx;
        int k_prev = (m_nx - 2) + j * m_nx;
        
        m_u[k_out] = m_u[k_prev];              
        m_v[k_out] = m_v[k_prev];              
    }
}

void Stokes::applyBoundaryToStar(float U_in) {
    // Bottom (j=0) and Top (j=ny-1): Solid walls (No-slip)
    for (int i = 0; i < m_nx; ++i) {
        m_uStar[i + 0 * m_nx] = 0.0;            // Bottom
        m_vStar[i + 0 * m_nx] = 0.0;           
        m_uStar[i + (m_ny - 1) * m_nx] = 0.0;   // Top
        m_vStar[i + (m_ny - 1) * m_nx] = 0.0;
    }

// Left (i=0): Inlet (Fluid entering)
    for (int j = 0; j < m_ny; ++j) {
        m_uStar[0 + j * m_nx] = U_in;
        m_vStar[0 + j * m_nx] = 0.0;
    }
    
    // Right (i=nx-1): Outlet (Zero-gradient / free exit)
    for (int j = 0; j < m_ny; ++j) {
        int k_out = (m_nx - 1) + j * m_nx;
        int k_prev = (m_nx - 2) + j * m_nx;
        
        m_uStar[k_out] = m_uStar[k_prev];  // Copy velocity from the  
        m_vStar[k_out] = m_vStar[k_prev];  // cell to the left
    }
}

void Stokes::correct(const std::vector<float>& p, float dt, float rho, float dxi, float dyi) {
    
    // Cylinder parameters
    float dx = 1.0 / dxi;
    float dy = 1.0 / dyi;
    float cx = 0.3;      // X-center of cylinder (30% down the tunnel)
    float cy = 0.49f;      // Y-center of cylinder (Middle of tunnel)
    float R = 0.1;       // Radius of the cylinder
    float R2 = R * R;

    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) {
            
            // 1D index for 2D grid
            int k = i + j * m_nx;

            // 1. Correct the velocity using the pressure gradient
            m_u[k] = m_uStar[k] - (dt / rho) * (p[k+1] - p[k-1]) * 0.5 * dxi;
            m_v[k] = m_vStar[k] - (dt / rho) * (p[k+m_nx] - p[k-m_nx]) * 0.5 * dyi;

            // 2. Calculate physical (x, y) coordinates of this cell
            // (Matching your Mesh.cpp formulation: (i - 0.5) * dx )
            float x = (i - 0.5) * dx;
            float y = (j - 0.5) * dy;

            // 3. CYLINDER MASK: If cell is inside the circle, force velocity to zero
            if ( ((x - cx) * (x - cx) + (y - cy) * (y - cy)) <= R2 ) {
                m_u[k] = 0.0;
                m_v[k] = 0.0;
            }
        }
    }
}

void Stokes::exportFrame(string filename, int frame, const vector<float>& p ) {

    ofstream out(filename, ios::binary);

    out.write((char*)m_u.data(), m_u.size() * sizeof(float));
    out.write((char*)m_v.data(), m_v.size() * sizeof(float));
    out.write((char*)p.data(), p.size() * sizeof(float));
    
    out.close();
}
