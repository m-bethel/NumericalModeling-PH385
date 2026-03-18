/*
Main File

 Project_9 main file: finding the enrgyenergy levels of quantum anharmonic
 oscillators.
 
 A note on units: we are assuming hbar=1 and m=1. This results in length
 units that are really funky and, more importantly, energy units that are
 essentially in units of hbar*omega (convenient for comparing n=2 levels
 to the analytic solutions for the harmonic oscillator).

(README)

Author: Miles Bethel (miles.d.bethel@gmail.com)
Date: 03/16/2026
*/
#include <iostream>
#include <limits>
#include "Solve.h"

using namespace std;

int main()
{
    int nmax = 10;
    int lmax = 10;
    double dx = 0.005;
    double dE = 0.02;
    double tol = 1E-6;
    double bdiv = 10.0;

    Solve solvEnergy (nmax, lmax, dx, dE, tol, bdiv);

    solvEnergy.output();

    return 0;
}
