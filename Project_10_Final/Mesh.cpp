/*
Mesh.cpp
-Sets up the 2D mesh for CFD flow

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#include "Mesh.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

Mesh::Mesh(int nx, int ny, double Lx, double Ly, int imin, int imax, int jmin, int jmax,
           double dx, double dy, double dxi, double dyi)

    : m_nx(nx), m_ny(ny),
      m_Lx(Lx), m_Ly(Ly),
      m_imin(imin), m_imax(imax),
      m_jmin(jmin), m_jmax(jmax),
      m_dx(dx), m_dy(dy),
      m_dxi(dxi), m_dyi(dyi){}

void Mesh::setUpMesh()
{
    //Allocate arrays with some extra space for "ghost cells"
    m_x.resize(m_imax + 2); // x(imin : imax+1) allocate a bit more for safety
    m_y.resize(m_jmax + 2);
    m_xm.resize(m_imax + 1);
    m_ym.resize(m_jmax + 1);

    // Create Uniform field lines
    double dx_step = m_Lx / m_nx;
    for (int i = m_imin; i <= m_imax +1; i++){
        m_x[i] = (i - m_imin) * dx_step;
    }
    double dy_step = m_Ly / m_ny;
    for (int j = m_jmin; j <= m_jmax +1; j++){
        m_y[j] = (j - m_jmin) * dy_step;
    }

    // Cell Center
    for (int i = m_imin; i <= m_imax; i++){
        m_xm[i] = 0.5 * (m_x[i] + m_x[i + 1]);
    }
    for (int j = m_jmin; j <= m_jmax; j++){
        m_ym[j] = 0.5 * (m_y[j] + m_y[j + 1]);
    }

    // Mesh Sizes
    m_dx = m_x[m_imin + 1] - m_x[m_imin];
    m_dy = m_y[m_jmin + 1] - m_y[m_jmin];

    m_dxi = 1.0 / m_dx;
    m_dyi = 1.0 / m_dy;

}

int Mesh::getNx() const
{
    return m_nx;
}
int Mesh::getNy() const
{
    return m_ny;
}
double Mesh::getDxi() const
{
    return m_dxi;
}
double Mesh::getDyi() const
{
    return m_dyi;
}