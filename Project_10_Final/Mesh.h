/*
Mesh.h 
- Sets up the initial Mesh

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#ifndef MESH_H
#define MESH_H

#include <string>
#include <vector>

class Mesh
{
    private:
        int m_nx, m_ny;
        int m_imin, m_imax, m_jmin, m_jmax;
        double m_Lx, m_Ly;
        double m_dx, m_dy, m_dxi, m_dyi;
        std::vector<double> m_x, m_y, m_xm, m_ym;

    public:

        Mesh(int nx, int ny, double Lx, double Ly, int imin, int imax, int jmin, int jmax,
             double dx, double dy, double dxi, double dyi);

        void setUpMesh();
        
        int getNx() const { return m_nx; }
        int getNy() const { return m_ny; }
        double getDxi() const { return m_dxi; }
        double getDyi() const { return m_dyi; }
};

#endif //MESH_H