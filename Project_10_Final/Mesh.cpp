/*
Mesh.cpp
-Sets up the 2D mesh for CFD flow

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#include "Mesh.h"

Mesh::Mesh(int nx, int ny, float Lx, float Ly)
    : m_nx(nx), m_ny(ny),
      m_Lx(Lx), m_Ly(Ly)
{
    // Calculate physical spacing
    m_dx = m_Lx / (float)(m_nx - 1);
    m_dy = m_Ly / (float)(m_ny - 1);
    
    // Pre-calculate inverse spacing for fast multiplication
    m_dxi = 1.0 / m_dx;
    m_dyi = 1.0 / m_dy;
}