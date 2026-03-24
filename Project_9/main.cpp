/*
Main File

 Project_9 main file: finding the energy levels of quantum anharmonic
 oscillators.
 
 A note on units: we are assuming hbar=1 and m=1. This results in length
 units that are really funky and, more importantly, energy units that are
 essentially in units of hbar*omega (convenient for comparing n=2 levels
 to the analytic solutions for the harmonic oscillator).

(README)



Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/20/2026
*/
#include <iostream>
#include <limits>
#include "Solve.h"

using namespace std;

int main()
{
    int nmax = 10;
    int lmax = 11;
    double dx = 0.005;
    double dE = 0.1;
    double tol = 1E-6;
    double bdiv = 10.0;
    double length = 10.0;
    

    Solve solver (lmax, dx, dE, tol, bdiv, length); //

    for (int n = 2; n <= nmax; n+=2)
    {
        double E = 0.1;
        for ( int state=0; state < lmax; state++){
            bool even_parity = (state % 2 == 0);

            E = solver.shooting(E, even_parity, n);

            string filename = "shoot_n" + to_string(n) + ".csv";
            solver.output(filename, state, n);
            E += 0.5;
        }
    }
    
    for (int n = 2; n <= nmax; n+=2)
    {
        double E = 0.1;
        for ( int state=0; state < lmax; state++){
            bool even_parity = (state % 2 == 0);

            E = solver.matching(E, even_parity, n);

            string filename = "match_n" + to_string(n) + ".csv";
            solver.output(filename, state, n);
            E += 0.5;
        }
    }
    return 0;
}
