/*
Mesh.h 
- Sets up the initial Mesh

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#ifndef MESH_H
#define MESH_H

class Mesh {
private:
    int m_nx, m_ny;
    float m_Lx, m_Ly;
    float m_dx, m_dy;   // Physical spacing
    float m_dxi, m_dyi; // Inverse spacing

public:
    Mesh(int nx, int ny, float Lx, float Ly);
    
    // Grid dimensions
    int getNx() const { return m_nx; }
    int getNy() const { return m_ny; }

    // Domain sizes
    float getLx() const { return m_Lx; }
    float getLy() const { return m_Ly; }

    // Physical spacing (Used for placing the cylinder)
    float getDx() const { return m_dx; }
    float getDy() const { return m_dy; }

    // Inverse spacing (Used for fast math in Poisson/Stokes)
    float getDxi() const { return m_dxi; }
    float getDyi() const { return m_dyi; }
};

#endif //MESH_H