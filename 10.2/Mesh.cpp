/*
Mesh.cpp
-Sets up the 2D mesh for CFD flow

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#include "Mesh.h"
#include <iostream>
#include <vector>

using namespace std;

Mesh::Mesh(int imin, int jmin, int imax, int jmax, 
             double Lx, double Ly, int nx, int ny)
           
    : m_nx(nx), m_ny(ny),
      m_Lx(Lx), m_Ly(Ly),
      m_imin(imin), m_jmin(jmin),
      m_imax(imax), m_jmax(jmax)
{
    setUpMesh();
}

void Mesh::setUpMesh()
{
    x.assign(m_imax + 2, 0.0);
    y.assign(m_jmax + 2, 0.0);
    xm.assign(m_imax + 1, 0.0);
    ym.assign(m_jmax + 1, 0.0);

    int total_cells = (m_imax + 2) * (m_jmax + 2);
    u.assign(total_cells, 0.0);
    v.assign(total_cells, 0.0);
    p.assign(total_cells, 0.0);

    double dx_spacing = m_Lx / m_nx;
    for (int i = 0; i <= m_ny; ++i){
        x[m_imin + i] = i * dx_spacing;
    }
    double dy_spacing = m_Ly / m_ny;
    for (int j = 0; j <= m_ny; ++j){
        y[m_imin + j] = j * dy_spacing;
    }
    for (int i = m_imin; i <= m_imax; ++i) {
        xm[i] = 0.5 * (x[i] + x[i + 1]);
    }
    for (int j = m_jmin; j <= m_jmax; ++j) {
        ym[j] = 0.5 * (y[j] + y[j + 1]);
    }
    double dx = x[m_imin + 1] - x[m_imin];
    double dy = y[m_jmin + 1] - y[m_jmin];
    
    double dxi = 1.0 / dx;
    double dyi = 1.0 / dy;
}

