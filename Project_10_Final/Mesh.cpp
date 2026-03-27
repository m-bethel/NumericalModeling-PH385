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

Mesh::Mesh(int nx, int ny, double Lx, double Ly)
           
    : m_nx(nx), m_ny(ny),
      m_Lx(Lx), m_Ly(Ly)
      {
        m_dx = m_Lx / (double) m_nx;
        m_dy = m_Ly / (double) m_ny;
        m_dxi = 1.0 / m_dx;
        m_dyi = 1.0 / m_dy;
      }

void Mesh::setUpMesh()
{
    //Allocate arrays with some extra space for "ghost cells"
    m_x.resize(m_nx);
    m_y.resize(m_ny);

    for (int i = 0; i < m_nx; ++i) {
        m_x[i] = (i - 0.5) * m_dx;
    }
    for (int j = 0; j < m_ny; ++j) {
        m_y[j] = (j - 0.5) * m_dy;
    }
}
