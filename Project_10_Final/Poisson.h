/*
Poisson.h 
- Defines the interface for the Poisson solver.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#ifndef POISSON_H
#define POISSON_H

#include "Mesh.h"
#include <vector>

class Poisson
{
    private:

        int m_nx, m_ny;
        int m_size;
        std::vector<double> m_values;   // Non-Zero entries
        std::vector<int> m_colIndices;  // Column index
        std::vector<int> m_rowPtr;      // Start Index
        std::vector<double> m_diag;     // for the fast solver
        std::vector<double> m_p;        // Pressure solution vector
        std::vector<double> m_RHS;      // Right-Hand Side vector (R)

    public:
    
        Poisson(const Mesh& mesh); // Constructor to build L
        void solve();
        std::vector<double> getP() const { return m_p; }

        void buildRHS(const std::vector<double>& uStar, 
                      const std::vector<double>& vStar, 
                      double dt, double rho, double dxi, double dyi);

};

#endif //POISSON_H