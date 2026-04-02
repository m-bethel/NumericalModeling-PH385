/*
Mesh.cpp
-Sets up the 2D mesh for CFD flow

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#include "Mesh.h"

// Constructor: Automatically called when a Mesh object is created in main.cpp
Mesh::Mesh(int nx, int ny, float Lx, float Ly)
    : m_nx(nx), m_ny(ny), m_Lx(Lx), m_Ly(Ly) // Initialize member variables
{
    // Calculate the physical size of one grid cell. 
    // We subtract 1 because N points create N-1 spaces between them.
    m_dx = m_Lx / (float)(m_nx - 1);
    m_dy = m_Ly / (float)(m_ny - 1);
    
    // Division operations are slow for a CPU. 
    // By storing (1/dx) here, we can use fast multiplication everywhere else in the code.
    m_dxi = 1.0f / m_dx;
    m_dyi = 1.0f / m_dy;
}