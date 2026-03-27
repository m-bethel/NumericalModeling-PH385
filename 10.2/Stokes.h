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

    double m_nu; // Viscosity
    int m_nx, m_ny; // grid parameter
    int imin, jmin, imax, jmax;

    std::vector<double> m_u, m_v;       // Current velocities (n)
    std::vector<double> m_uStar, m_vStar; // Intermediate velocities (*)

public:
    Stokes(int nx, int ny, const Mesh& mesh, double nu);
    
    // The Predictor Step (Sections 4 & 5)
    void predict(double dt, double dxi, double dyi);
    
    // The Corrector Step (Section 7)
    void correct(const std::vector<double>& p, double dt, double rho, double dxi, double dyi);

    // Boundary Conditions
    void applyBoundary(double U_lid);
    void applyBoundaryToStar(double U_lid);

    void exportFrame(std::string filename, int frame, const std::vector<double>& p);

    // Getters for the Poisson RHS
    const std::vector<double>& getUStar() const { return m_uStar; }
    const std::vector<double>& getVStar() const { return m_vStar; }
};

#endif