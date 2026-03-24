/*
Solve.cpp


Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/

#include "Solve.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

Solve::Solve( int lmax, double dx, double dE, 
              double tol, double bdiv, double length)
    : m_lmax(lmax),
      m_dx(dx),
      m_dE(dE),
      m_tol(tol),
      m_bdiv(bdiv),
      m_length(length){}

// Calculates the potential well shape based on the power 'n'
double Solve::calculatePotential(double x, int n_power)
{
    double result = 1.0;
    double x_squared = x*x;
    for (int i = 0; i < n_power /2 ; i++) result *= x_squared;
    return result/n_power;
}

// THE SHOOTING METHOD: Integrates from x=0 toward x=L
double Solve::shooting(double E, bool even_parity, int n_power)
{
    double deltaE = m_dE;
    double last_sign = 0.0;
    int steps = m_length / m_dx;

    while (abs(deltaE) > m_tol) {  
        m_psi.clear();
        // Set initial conditions at x=0 based on parity
        if(even_parity) {
            m_psi.push_back(1.00); // psi(0) = 1 (even symmetry)
            m_psi.push_back(1.00); // psi(dx) approx psi(0)
        } else {
            m_psi.push_back(-m_dx); // Slope at origin for odd states
            m_psi.push_back(0.00);  // psi(0) = 0 (odd symmetry)
        }

        for (int i = 1; i < steps - 1; ++i) {
            double x = (i - 1) * m_dx;
            double V_n = calculatePotential(x, n_power);
            // Numerov recurrence: calculates next psi based on previous two
            double next_psi = 2*m_psi[i] - m_psi[i-1] - 2*(m_dx*m_dx)*(E - V_n)*m_psi[i];
            m_psi.push_back(next_psi);

            // Stop if the wavefunction blows up (divergence)
            if (abs(m_psi.back()) > m_bdiv) break;
        }

        // Adjust energy guess based on whether the tail crossed the axis
        double current_sign = (m_psi.back() > 0) ? 1.0 : -1.0;
        if (last_sign != 0.00 && current_sign != last_sign) {
            deltaE = -deltaE / 2.0; // Binary search refinement
        }
        E += deltaE;
        last_sign = current_sign;
    }
    return E;
}

// THE MATCHING METHOD: Integrates from a stable boundary (L) back to the origin (0)
double Solve::matching(double E, bool even_parity, int n_power)
{
    double deltaE = m_dE;
    double last_sign = 0.0;
    
    // Higher power potentials are narrower; we shorten the integration range for stability
    double effective_length = (n_power >= 8) ? 4.0 : ((n_power > 4) ? 5.0 : m_length);
    int steps = effective_length / m_dx;

    while (abs(deltaE) > m_tol) {
        m_psi.clear();
        m_psi.push_back(0.0);    // Seed: psi(L) = 0 (must be zero at infinity)
        m_psi.push_back(1e-15);  // Small non-zero value to start the integration

        for (int i = 1; i < steps; i++) {
            double x = effective_length - (i * m_dx); // Moving backward toward 0
            double V_n = calculatePotential(x, n_power);
            double next_psi = 2 * m_psi[i] - m_psi[i-1] - 2 * (m_dx * m_dx) * (E - V_n) * m_psi[i];
            m_psi.push_back(next_psi);
        }

        // Check boundary conditions at the origin (x=0)
        double current_sign;
        if (even_parity) {
            // Even states must have zero slope at x=0
            current_sign = (m_psi.back() - m_psi[m_psi.size()-2] > 0) ? 1.0 : -1.0;
        } else {
            // Odd states must cross exactly at 0
            current_sign = (m_psi.back() > 0) ? 1.0 : -1.0;
        }

        if (last_sign != 0.0 && current_sign != last_sign) {
            deltaE = -deltaE / 2.0;
        }
        E += deltaE;
        last_sign = current_sign;
    }
    return E;
}

// Handles outputting results in a format Python can easily mirror
void Solve::output(string filename, int state, int n_power)
{
    ofstream outFile(filename, std::ios::app);
    if (!outFile.is_open()) return;

    // Detection logic to see if this was a backward (matching) run
    bool is_matching = (!m_psi.empty() && m_psi[0] == 0.0);
    double effective_length = (n_power >= 8) ? 4.0 : ((n_power > 4) ? 5.0 : m_length);

    if (is_matching) {
        // To keep x-values ascending for the plotting script, we write matching data in reverse
        for (int i = m_psi.size() - 1; i >= 0; i--) {
            double dist_from_origin = (m_psi.size() - 1 - i) * m_dx;
            outFile << state << "," << dist_from_origin << "," << m_psi[i] << "\n";
        }
    } else {
        for (int i = 0; i < m_psi.size(); i++) {
            double x = i * m_dx;
            outFile << state << "," << x << "," << m_psi[i] << "\n";
        }
    }
    outFile.close();
}