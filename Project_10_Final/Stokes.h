/*
Stokes.h & Stokes.cpp
- Solves the Navier-Stokes equations to figure out how the fluid's velocity 
    changes over time.
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
    std::vector<float> m_u, m_v;         // The actual, physical horizontal (u) and vertical (v) velocities
    std::vector<float> m_uStar, m_vStar; // The "predicted" velocities before pressure is applied (*)
    float m_nu;                          // Viscosity of the fluid
    int m_nx, m_ny;                      // Grid dimensions copied from Mesh

public:
    Stokes(const Mesh& mesh, float nu);
    
    // Calculates how momentum carries the fluid forward
    void predict(float dt, float dxi, float dyi);
    
    // Uses the pressure field to enforce conservation of mass
    void correct(const std::vector<float>& p, float dt, float rho, float dxi, float dyi);

    // Forces velocities at the edges of the simulation to behave realistically
    void applyBoundary(float U_in);
    void applyBoundaryToStar(float U_in);

    // Dumps the arrays to a binary file for Python to read
    void exportFrame(std::string filename, int frame, const std::vector<float>& p);

    // Allows the Poisson class to read the predicted velocities
    const std::vector<float>& getUStar() const { return m_uStar; }
    const std::vector<float>& getVStar() const { return m_vStar; }
};

#endif