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

// Constructor: Sets up the memory arrays for the velocity grids
Stokes::Stokes(const Mesh& mesh, float nu) : m_nu(nu) {
    m_nx = mesh.getNx();
    m_ny = mesh.getNy();
    int size = m_nx * m_ny; // Total number of cells

    // Allocate memory and fill with 0.0 (Starting with perfectly still fluid)
    m_u.assign(size, 0.0f);
    m_uStar.assign(size, 0.0f);
    m_v.assign(size, 0.0f);
    m_vStar.assign(size, 0.0f);
}

void Stokes::predict(float dt, float dxi, float dyi) {
    // Pre-calculate squared values for diffusion math
    float dxi2 = dxi * dxi;
    float dyi2 = dyi * dyi;
    float dx = 1.0f / dxi;
    float dy = 1.0f / dyi;

    // --- CYLINDER OBSTACLE DEFINITION ---
    // We place a circle slightly off-center (0.49 instead of 0.5) to break 
    // the mathematical symmetry and trigger the vortex street.
    float cx = 0.3f;     // Center X (30% down the tunnel)
    float cy = 0.49f;    // Center Y (Slightly off perfect center)
    float R = 0.1f;      // Radius of the cylinder
    float R2 = R * R;    // Radius squared (for faster distance math)

    // --- U-VELOCITY (HORIZONTAL) PREDICTION ---
    // Loop over every cell in the domain (ignoring the absolute outer edges)
    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) { 
            int k = i + j * m_nx; // Convert 2D (x,y) to 1D memory array index

            // Calculate the physical (x, y) coordinate of this specific cell
            float x = i * dx; 
            float y = (j + 0.5f) * dy;

            // Check if this cell is inside the solid cylinder using Pythagorean theorem (x^2 + y^2 = R^2)
            if (((x - cx) * (x - cx) + (y - cy) * (y - cy)) <= R2) {
                m_uStar[k] = 0.0f; // Inside a solid object, fluid cannot move
            } else {
                // ADVECTION: How the fluid's own velocity carries it forward.
                // We use an "Upwind" scheme: We look at the fluid coming *from behind* us.
                float v_here = 0.25f * (m_v[k] + m_v[k + 1] + m_v[k - m_nx] + m_v[k + 1 - m_nx]); // Average vertical speed
                float u_grad_x = (m_u[k] > 0) ? (m_u[k] - m_u[k-1]) * dxi : (m_u[k+1] - m_u[k]) * dxi; // X-gradient
                float u_grad_y = (v_here > 0) ? (m_u[k] - m_u[k-m_nx]) * dyi : (m_u[k+m_nx] - m_u[k]) * dyi; // Y-gradient
                float convection = m_u[k] * u_grad_x + v_here * u_grad_y;
                
                // DIFFUSION: How viscosity causes fluid to pull on neighboring fluid (like honey).
                // Calculated using the 2nd derivative (Laplacian) of velocity.
                float diffusion = m_nu * (
                    (m_u[k+1] - 2.0f*m_u[k] + m_u[k-1]) * dxi2 +
                    (m_u[k+m_nx] - 2.0f*m_u[k] + m_u[k-m_nx]) * dyi2
                );

                // Calculate the final predicted velocity
                m_uStar[k] = m_u[k] + dt * (-convection + diffusion);
            }
        }
    }

    // --- V-VELOCITY (VERTICAL) PREDICTION ---
    // (This follows the exact same logical steps as the U-velocity, but vertically)
    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) {
            int k = i + j * m_nx;
            float x = (i + 0.5f) * dx;
            float y = j * dy;

            if (((x - cx) * (x - cx) + (y - cy) * (y - cy)) <= R2) {
                m_vStar[k] = 0.0f;
            } else {
                float u_here = 0.25f * (m_u[k] + m_u[k - 1] + m_u[k + m_nx] + m_u[k + m_nx - 1]);
                float v_here = m_v[k];
                float v_grad_x = (u_here > 0) ? (m_v[k] - m_v[k-1]) * dxi : (m_v[k+1] - m_v[k]) * dxi;
                float v_grad_y = (m_v[k] > 0) ? (m_v[k] - m_v[k-m_nx]) * dyi : (m_v[k+m_nx] - m_v[k]) * dyi;
                
                float convection = v_here * v_grad_y + u_here * v_grad_x;
                float diffusion = m_nu * (
                    (m_v[k+1] - 2.0f*m_v[k] + m_v[k-1]) * dxi2 +
                    (m_v[k+m_nx] - 2.0f*m_v[k] + m_v[k-m_nx]) * dyi2
                );

                m_vStar[k] = m_v[k] + dt * (-convection + diffusion);
            }
        }
    }
}

// Applies physics rules to the outer edges of the simulation box
void Stokes::applyBoundary(float U_in){
    // Bottom and Top Walls: "No-Slip Condition"
    // Friction forces the fluid touching the wall to have zero velocity.
    for (int i = 0; i < m_nx; ++i) {
        m_u[i + 0 * m_nx] = 0.0f;               
        m_v[i + 0 * m_nx] = 0.0f;           
        m_u[i + (m_ny - 1) * m_nx] = 0.0f;       
        m_v[i + (m_ny - 1) * m_nx] = 0.0f;
    }

    // Left Wall: "Inlet Condition"
    // We force fluid into the tunnel at a constant speed (U_in).
    for (int j = 0; j < m_ny; ++j) {
        m_u[0 + j * m_nx] = U_in;               
        m_v[0 + j * m_nx] = 0.0f;
    }
    
    // Right Wall: "Outlet Condition"
    // We let fluid flow out freely by copying the velocity from the cell right next to it.
    for (int j = 0; j < m_ny; ++j) {
        int k_out = (m_nx - 1) + j * m_nx;
        int k_prev = (m_nx - 2) + j * m_nx;
        m_u[k_out] = m_u[k_prev];              
        m_v[k_out] = m_v[k_prev];              
    }
}

// Exactly the same as above, but applied to the temporary "Predicted" arrays
void Stokes::applyBoundaryToStar(float U_in) {
    for (int i = 0; i < m_nx; ++i) {
        m_uStar[i + 0 * m_nx] = 0.0f;            
        m_vStar[i + 0 * m_nx] = 0.0f;           
        m_uStar[i + (m_ny - 1) * m_nx] = 0.0f;   
        m_vStar[i + (m_ny - 1) * m_nx] = 0.0f;
    }
    for (int j = 0; j < m_ny; ++j) {
        m_uStar[0 + j * m_nx] = U_in;
        m_vStar[0 + j * m_nx] = 0.0f;
    }
    for (int j = 0; j < m_ny; ++j) {
        int k_out = (m_nx - 1) + j * m_nx;
        int k_prev = (m_nx - 2) + j * m_nx;
        m_uStar[k_out] = m_uStar[k_prev];  
        m_vStar[k_out] = m_vStar[k_prev];  
    }
}

void Stokes::correct(const std::vector<float>& p, float dt, float rho, float dxi, float dyi) {
    
    // Cylinder parameters (Must exactly match the predictor step)
    float dx = 1.0f / dxi;
    float dy = 1.0f / dyi;
    float cx = 0.3f;      
    float cy = 0.49f;      
    float R = 0.1f;       
    float R2 = R * R;

    for (int j = 1; j < m_ny - 1; ++j) {
        for (int i = 1; i < m_nx - 1; ++i) {
            int k = i + j * m_nx;

            // CORRECTOR MATH: Fluid naturally flows from High Pressure to Low Pressure.
            // We calculate the pressure difference (gradient) across the cell, and subtract 
            // that from our predicted velocity to get the final, true velocity.
            m_u[k] = m_uStar[k] - (dt / rho) * (p[k+1] - p[k-1]) * 0.5f * dxi;
            m_v[k] = m_vStar[k] - (dt / rho) * (p[k+m_nx] - p[k-m_nx]) * 0.5f * dyi;

            float x = (i - 0.5f) * dx;
            float y = (j - 0.5f) * dy;

            // Apply the cylinder mask one final time to ensure no fluid leaked inside
            if ( ((x - cx) * (x - cx) + (y - cy) * (y - cy)) <= R2 ) {
                m_u[k] = 0.0f;
                m_v[k] = 0.0f;
            }
        }
    }
}

// Saves the entire memory block of data to a pure binary file for fast Python rendering
void Stokes::exportFrame(string filename, int frame, const vector<float>& p ) {
    ofstream out(filename, ios::binary);
    // Writes raw bytes directly from RAM to the Hard Drive
    out.write((char*)m_u.data(), m_u.size() * sizeof(float));
    out.write((char*)m_v.data(), m_v.size() * sizeof(float));
    out.write((char*)p.data(), p.size() * sizeof(float));
    out.close();
}