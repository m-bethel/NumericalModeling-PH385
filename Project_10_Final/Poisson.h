/*
Poisson.h 
- Defines the interface for the Poisson solver.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/


#ifndef POISSON_H
#define POISSON_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <omp.h>
#include "Mesh.h"

class Poisson {
public:
    Poisson(const Mesh& mesh);
    
    void solve(int max_Iterations, float tol, float dxi, float dyi);
    void buildRHS(const std::vector<float>& uStar, const std::vector<float>& vStar, 
                  float dt, float rho, float dxi, float dyi);
    
    const std::vector<float>& getP() const { return m_p; }

private:
    int m_nx, m_ny;
    float m_dx, m_dy;
    std::vector<float> m_p;
    std::vector<float> m_RHS;
    std::vector<bool> m_isSolid; // The "Fast" Mask
};

#endif //POISSON_H