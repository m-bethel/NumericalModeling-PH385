/*
Solve.h 
- Defines the interface for the Schrodinger equation solver.

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#ifndef SOLVE_H
#define SOLVE_H

#include <string>
#include <vector>

class Solve
{
    private:
        int m_lmax;          // Maximum number of energy states to find
        double m_dx;        // Step size for the spatial grid (x)
        double m_dE;        // Initial energy step for the shooting/matching method
        double m_tol;       // Convergence tolerance for finding energy eigenvalues
        double m_bdiv;      // Divergence threshold for the shooting method
        double m_length;    // Maximum spatial extent (L) for integration
        std::vector<double> m_psi; // Storage for the calculated wavefunction values
    
    public:
        // Constructor to initialize simulation parameters
        Solve( int lmax, double dx, double dE,
               double tol, double bdiv, double length);

        // Finds energy levels by integrating from a boundary back to the origin
        double matching(double E, bool even_parity, int n_power);
        
        // Finds energy levels by integrating from the origin out to a boundary
        double shooting(double E, bool even_parity, int n_power);
        
        // Calculates V(x) = x^n / n potential energy at a given point
        double calculatePotential(double x, int n_power);
        
        // Writes state data [state, x, psi] to a CSV file for plotting
        void output(std::string filename, int state, int n_power);

};

#endif //SOLVE_H