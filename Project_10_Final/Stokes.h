/*
Stokes.h 
- Defines the interface for the Stokes solver.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#ifndef STOKES_H
#define STOKES_H

#include <vector>
#include <string>
#include "Mesh.h"

class Stokes {
private:
    std::vector<float> m_u, m_v;       // Current velocities (n)
    std::vector<float> m_uStar, m_vStar; // Intermediate velocities (*)
    float m_nu; // Viscosity
    int m_nx, m_ny;

public:
    Stokes(const Mesh& mesh, float nu);
    
    // The Predictor Step (Sections 4 & 5)
    void predict(float dt, float dxi, float dyi);
    
    // The Corrector Step (Section 7)
    void correct(const std::vector<float>& p, float dt, float rho, float dxi, float dyi);

    // Boundary Conditions
    void applyBoundary(float U_in);
    void applyBoundaryToStar(float U_in);

    void exportFrame(std::string filename, int frame, const std::vector<float>& p);

    // Getters for the Poisson RHS
    const std::vector<float>& getUStar() const { return m_uStar; }
    const std::vector<float>& getVStar() const { return m_vStar; }
};

#endif