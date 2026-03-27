/*
Mesh.h 
- Sets up the initial Mesh

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#pragma once
#include <vector>

class Mesh
{
    private:
    //grid parameters
        int m_nx, m_ny, m_imin, m_jmin, m_imax, m_jmax;
        double m_Lx, m_Ly;
    //Mesh sizing variables
        double dx, dy, dxi, dyi;
    //Mesh Vectors
        std::vector<double> x, y;
        std::vector<double> xm, ym;
        std::vector<double> u; // x-velocity
        std::vector<double> v; // y-velocity
        std::vector<double> p; // pressure

    public:

        Mesh(int imin, int jmin, int imax, int jmax, 
             double Lx, double Ly, int nx, int ny);

        void setUpMesh();\

        //Utilizing In-lining
        int ID(int i, int j) const {
            return i + j *(m_imax + 2);
        }
        
        //int getNx() const { return m_nx; }
        //int getNy() const { return m_ny; }
        //double getDxi() const { return m_dxi; }
        //double getDyi() const { return m_dyi; }
        //double getX(int i) const { return m_x[i]; }
        //double getY(int j) const { return m_y[j]; }
};