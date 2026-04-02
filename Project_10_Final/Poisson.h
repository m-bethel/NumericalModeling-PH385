/*
Poisson.h & Poisson.cpp
- Solves the Pressure Poisson Equation. This forces the fluid to obey Conservation 
  of Mass (making it incompressible) by calculating the pressure field needed to 
  push overlapping fluid apart.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#ifndef POISSON_H
#define POISSON_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <omp.h> // OpenMP for Multi-Core CPU Parallel Processing
#include "Mesh.h"

class Poisson {
public:
    Poisson(const Mesh& mesh);
    
    // The iterative math solver
    void solve(int max_Iterations, float tol, float dxi, float dyi);
    
    // Calculates where the fluid is bunching up based on the Predictor step
    void buildRHS(const std::vector<float>& uStar, const std::vector<float>& vStar, 
                  float dt, float rho, float dxi, float dyi);
    
    const std::vector<float>& getP() const { return m_p; }

private:
    int m_nx, m_ny;
    float m_dx, m_dy;
    std::vector<float> m_p;       // The Pressure field
    std::vector<float> m_RHS;     // The "Right Hand Side" (Divergence / bunching up)
    std::vector<bool> m_isSolid;  // A pre-calculated mask of the cylinder to speed up math
};

#endif //POISSON_H