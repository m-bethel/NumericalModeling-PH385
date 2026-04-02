/*
Mesh.h & Mesh.cpp
- Defines the physical and computational boundaries of our simulation.
- Sets up the initial Mesh

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/


#ifndef MESH_H
#define MESH_H

class Mesh {
private:
    int m_nx, m_ny;       // Number of cells (grid resolution)
    float m_Lx, m_Ly;     // Physical dimensions of the simulation box
    float m_dx, m_dy;     // Physical width and height of a single cell
    float m_dxi, m_dyi;   // "Inverse" spacing (1/dx). We pre-calculate this because multiplication is faster than division for CPUs.

public:
    Mesh(int nx, int ny, float Lx, float Ly);
    
    // Standard Getters to allow other classes to read these values safely
    int getNx() const { return m_nx; }
    int getNy() const { return m_ny; }
    float getLx() const { return m_Lx; }
    float getLy() const { return m_Ly; }
    float getDx() const { return m_dx; }
    float getDy() const { return m_dy; }
    float getDxi() const { return m_dxi; }
    float getDyi() const { return m_dyi; }
};

#endif //MESH_H